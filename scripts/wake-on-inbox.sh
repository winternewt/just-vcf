#!/usr/bin/env bash
# Claude Code only: runs watch-inbox.sh detached and wakes the session through its inbox socket.
#
# Claude Code stops every background task it owns after a deadline (2.1.285: 30 min by default,
# 2 h at most; its Monitor tool: 30 min), so the one-shot arming in TRIAGE_LOOP.md § 2 dies on a
# timer there. This runs outside the harness and posts each settle as a user message to the
# session's own inbox socket, $CLAUDE_CODE_MESSAGING_SOCKET, which Claude Code exports to every hook
# and Bash command. An idle session starts a new turn when one arrives.
#
# The socket and its auth line are documented
# (https://code.claude.com/docs/en/cross-session-messaging); the message line is not. The shape used
# here, {"type":"user","message":{"role":"user","content":...}}, is taken from Claude Code's own help
# text. ARMED_NOTICE=1 posts one "armed" line first, and that line is the live check: if an update
# changes the format, the arming turn sees no notice, rather than a reporter's item waiting unseen.
#
# It keeps watching after each event, so nothing needs re-arming. It ends when the watcher ends
# (superseded by a newer arming, exit 3, or stopped), or when the socket is gone, which is what a
# closed session leaves behind. A /clear keeps the process and the socket, so it survives that.
# Every env knob of watch-inbox.sh (FILE, COOLDOWN, POLL, BRANCH, ...) passes through.
#
# Arm it from a hook (arm-on-prompt.py), never from the agent: auto mode's classifier refuses to let
# a session launch something that wakes itself, and a hook is not classified.
set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOCK=${CLAUDE_CODE_MESSAGING_SOCKET:?run from inside a Claude Code session}
TOKEN=${CLAUDE_CODE_MESSAGING_TOKEN:-}
READ_TIMEOUT=${READ_TIMEOUT:-60}   # how often a quiet loop checks the watcher and socket are alive

post() {
    SOCK=$SOCK TOKEN=$TOKEN TEXT=$1 python3 -c '
import json, os, socket
s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
s.connect(os.environ["SOCK"])
lines = [{"type": "auth", "token": os.environ["TOKEN"]}] if os.environ["TOKEN"] else []
lines.append({"type": "user", "message": {"role": "user", "content": os.environ["TEXT"]}})
s.sendall("".join(json.dumps(x) + "\n" for x in lines).encode())
s.close()
'
}

if [ "${ARMED_NOTICE:-}" = 1 ]; then
    post "[watcher $(date -u +%FT%TZ)] armed on ${FILE:-the inbox}" || exit 1
fi

coproc W { exec "$HERE/watch-inbox.sh"; }
p=$W_PID
# Keep our own copy of the read end: bash unsets W once the coprocess exits. A blocking read was
# seen to outlive a killed watcher, so it wakes every READ_TIMEOUT seconds to check both ends.
exec {rfd}<&"${W[0]}"
while :; do
    # Capture read's own status. `if ! read …; then [ $? … ]` sees the negation's 0, and turned
    # every quiet minute into a dead watcher.
    IFS= read -r -t "$READ_TIMEOUT" line <&"$rfd"
    rc=$?                           # >128 is the timeout; anything else non-zero is EOF
    if [ "$rc" -ne 0 ]; then
        [ "$rc" -gt 128 ] && kill -0 "$p" 2>/dev/null && [ -S "$SOCK" ] && continue
        break
    fi
    case "$line" in *"nothing pending"*|*paus*|*resum*) continue;; esac
    [ -S "$SOCK" ] || break
    # The timestamp keeps two identical settles apart: the inbox drops a quick repeat.
    post "[watcher $(date -u +%FT%TZ)] $line" || break
done
kill "$p" 2>/dev/null
wait "$p"
