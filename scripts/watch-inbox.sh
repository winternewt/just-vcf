#!/usr/bin/env bash
# Debounced watcher for the document-as-inbox triage loop. See TRIAGE_LOOP.md.
#
# Emits one line on stdout when the watched file has stopped changing for COOLDOWN seconds.
# Consecutive saves inside the cooldown collapse into a single event, because each mtime bump restarts
# the timer. Meant to be driven by whatever turns a line of stdout into a notification — an agent
# harness, `while read`, a desktop notifier — and it works piped to anything.
#
# It never fires for a change that predates it: `last` is seeded from the current mtime and `dirty`
# starts clear, so an already-settled edit stays quiet. Run the ledger once at startup to pick up
# whatever is already pending.
#
# COOLDOWN is sized for an agent author, not a human one. When reports are written through an agent the
# pattern is a burst of edits — five in a minute is normal — with gaps wherever the agent stops to read
# or probe something. A one-minute timer fires in the middle of such a run and triages half a report.
# 150s clears the pauses that show up in practice. The only cost of waiting is latency; the cost of
# firing early is a reply to a half-written item.
#
# Nothing needs installing: `stat` polling is enough at this cadence, and it works where inotify-tools,
# entr, fswatch and python watchdog are all absent.
#
# If your triage pass is permitted to commit, it should only watch while the tree is on the branch
# that permit covers: a feature branch — or a detached HEAD — is somebody's own work, and unattended
# triage commits would land on top of whatever they are mid-way through. Off that branch the watcher
# idles at BRANCH_PAUSE instead of POLL, says so once, and says so again when it resumes. It does not
# touch `last` while paused, so an edit written during the pause is still picked up on the way back
# rather than lost. Set BRANCH to your own name for it, or to the empty string to switch the whole
# behaviour off; outside a git work tree there is no branch to speak of and it never pauses.
#
# One watcher per watched file, and the newest wins. Arming it again replaces the running one instead
# of adding a second: the older watcher belongs to an earlier run of the agent, so it notifies nobody who
# is listening, and two live watchers report every settle twice. The key is the watched file's resolved
# path, so a sibling repo's watcher is a different key and is never touched. A pidfile names the owner.
# On start the watcher records itself there and stops the previous owner, but only after checking that
# pid's command line, because a recycled pid could belong to anything. Every poll it checks the pidfile
# and exits if it no longer owns it. That exit is what settles two watchers started in the same second,
# and it also covers an old watcher whose kill was refused. No flock is used, because macOS has none.
#
# A replaced watcher exits with SUPERSEDED (3), whether it noticed at a poll or was sent TERM by its
# successor. It is the one exit that is neither an event nor a fault, and whoever armed it should be able
# to tell it apart from both without reading anything: a wrapper that ends in `wait` on this process
# reports 3 as the task's status. A TERM or INT while this watcher still owns the pidfile is somebody
# stopping it on purpose, and exits 0. Any other non-zero status is a real failure.
#
# This is the only one of the three that is really bash. The ledger is Python and is invoked through
# $PYTHON below rather than as a bare path, so neither its exec bit nor its shebang is load-bearing
# (TRIAGE_LOOP.md § 5 — the extension gotcha).
#
#   FILE=<path> COOLDOWN=<seconds> POLL=<seconds> BRANCH=<name|empty> BRANCH_PAUSE=<seconds> \
#     RUNBOOK=<path> CAP=<n> ./watch-inbox.sh
set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# Local to just-vcf: FILE and RUNBOOK derive from this script's location rather than $PWD.
FILE=${FILE:-${INBOX:-$HERE/../docs/FEEDBACK.md}}
PYTHON=${PYTHON:-python3}
LEDGER=${LEDGER:-$HERE/triage-state.py}
RUNBOOK=${RUNBOOK:-docs/CONSUMER_TRIAGE_LOOP.md}
COOLDOWN=${COOLDOWN:-150}
POLL=${POLL:-10}
BRANCH=${BRANCH-main}          # `-`, not `:-`: BRANCH= is an explicit "never pause", not an unset
BRANCH_PAUSE=${BRANCH_PAUSE:-900}
CAP=${CAP:-8}

# The tree to ask about is the one holding the watched file, not the one holding this script — the
# three scripts may live anywhere, including outside the repo they serve.
REPO=$(cd "$(dirname "$FILE")" 2>/dev/null && pwd || echo "$PWD")
# Answered once: a directory does not become a work tree while the watcher runs, and if it is not one
# now then `branch` means nothing here and the pause must never fire.
if [ -n "$BRANCH" ] && git -C "$REPO" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    watch_branch=1
else
    watch_branch=0
fi

# The singleton. Keyed on the resolved path of the watched file, not on this script.
FILE_ABS=$(cd "$(dirname "$FILE")" 2>/dev/null && echo "$PWD/${FILE##*/}" || echo "$FILE")
PIDFILE=${PIDFILE:-${XDG_RUNTIME_DIR:-${TMPDIR:-/tmp}}/watch-inbox-$(printf %s "$FILE_ABS" | cksum | cut -d' ' -f1).pid}
old=$(cat "$PIDFILE" 2>/dev/null || true)
echo $$ >"$PIDFILE"
if [ -n "$old" ] && [ "$old" != $$ ] && kill -0 "$old" 2>/dev/null &&
   ps -o args= -p "$old" 2>/dev/null | grep -q "${BASH_SOURCE[0]##*/}"; then
    kill "$old" 2>/dev/null && echo "replaced watcher pid $old on ${FILE##*/}" >&2
fi
owner() { [ "$(cat "$PIDFILE" 2>/dev/null)" = $$ ]; }
SUPERSEDED=3
# Leave the pidfile behind only if it names somebody else. A replaced watcher must not delete its successor's.
trap 'owner && rm -f "$PIDFILE"' EXIT
# The successor writes the pidfile before it signals, so a TERM arriving here already reads as not-owner.
trap 'owner && exit 0; exit "$SUPERSEDED"' TERM INT
# Bash runs a trap only after its foreground child returns, so a plain `sleep 900` during a branch pause
# would keep a replaced watcher alive for fifteen minutes. `wait` is interruptible.
nap() { sleep "$1" & wait $!; }

mtime() { stat -c %Y "$FILE" 2>/dev/null || stat -f %m "$FILE" 2>/dev/null || echo 0; }
# A detached HEAD has no symbolic ref, and is no more a place to commit unattended than a branch is.
current_branch() { git -C "$REPO" symbolic-ref --quiet --short HEAD 2>/dev/null || echo "(detached HEAD)"; }

last=$(mtime)
dirty=0
paused=""

while true; do
    owner || exit "$SUPERSEDED"     # replaced by a newer arming
    if [ "$watch_branch" = 1 ]; then
        on=$(current_branch)
        if [ "$on" != "$BRANCH" ]; then
            # Announce the transition once. A pause nobody can see reads as a dead watcher, and
            # repeating it every quarter hour would be the noise the event filter exists to avoid.
            if [ "$paused" != "$on" ]; then
                paused=$on
                echo "${FILE##*/} watch paused: tree is on $on, not $BRANCH — a branch is the human's own work"
            fi
            nap "$BRANCH_PAUSE"
            continue
        fi
        if [ -n "$paused" ]; then
            paused=""
            echo "${FILE##*/} watch resumed: back on $BRANCH"
        fi
    fi

    nap "$POLL"
    owner || exit "$SUPERSEDED"
    now=$(mtime)

    if [ "$now" != "$last" ]; then
        last=$now
        dirty=1                     # still moving: restart the cooldown
        continue
    fi

    [ "$dirty" = 1 ] || continue
    [ $(($(date +%s) - last)) -ge "$COOLDOWN" ] || continue
    dirty=0

    # The event line is capped: with a 17-item backlog an uncapped one listed every single item.
    pending=$(INBOX="$FILE" "$PYTHON" "$LEDGER" "$FILE" --pending 2>/dev/null |
              awk -v cap="$CAP" 'NF { n++; if (n <= cap) { printf "%s%s(%s)", sep, $2, $1; sep = " " } }
                   END { if (n > cap) printf " +%d more", n - cap }')
    if [ -z "$pending" ]; then
        echo "${FILE##*/} settled: nothing pending"
    else
        # Naming the runbook matters: to an agent whose context was just cleared, this line is the
        # entire brief — it has to say where the algorithm is written down.
        echo "${FILE##*/} settled: $pending — triage per $RUNBOOK"
    fi
done
