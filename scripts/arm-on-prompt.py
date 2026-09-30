#!/usr/bin/env python3
"""Claude Code UserPromptSubmit hook: arm wake-on-inbox.sh when a prompt cites the runbook.

The triage seat is the session pointed at the runbook, so citing it in a prompt, in any wording, is
what arms the watcher, and the watcher belongs to that session and no other. RUNBOOK names the file to
look for (default TRIAGE_LOOP.md); FILE names the inbox, as it does for watch-inbox.sh. Register it in
the repo's committed .claude/settings.json:

    {"hooks": {"UserPromptSubmit": [{"hooks": [{"type": "command", "timeout": 10,
      "command": "FILE=\\"$CLAUDE_PROJECT_DIR/docs/INBOX.md\\" python3 \\"$CLAUDE_PROJECT_DIR/tools/arm-on-prompt.py\\""}]}]}}

Never blocks. Prints one line to stdout, which Claude Code adds to the turn, saying what it did.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

WRAPPER = Path(__file__).resolve().parent / "wake-on-inbox.sh"
WATCHER_POST = "[watcher "  # the prefix wake-on-inbox.sh puts on every line it posts


def main() -> int:
    try:
        data = json.load(sys.stdin)
    except ValueError:
        return 0
    prompt = str(data.get("prompt", ""))
    # The watcher's own posts cite the runbook, and whether Claude Code runs this hook on a socket
    # message is undocumented. A machine-written line re-arming the machine that wrote it is the
    # loop to rule out, not to bound, so none of its posts may arm anything.
    if WATCHER_POST in prompt or os.environ.get("RUNBOOK", "TRIAGE_LOOP.md") not in prompt:
        return 0
    sock = os.environ.get("CLAUDE_CODE_MESSAGING_SOCKET")
    if not sock:
        print("[triage watcher] NOT armed: this session exports no inbox socket. Use the one-shot arming.")
        return 0
    subprocess.Popen(
        [str(WRAPPER)],
        env={**os.environ, "ARMED_NOTICE": "1"},
        start_new_session=True,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    print(
        f"[triage watcher] armed on {sock}. A '[watcher …] armed' message should arrive this turn; "
        "if none does, the inbox format changed: say so, and use the one-shot arming."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
