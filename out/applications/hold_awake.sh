#!/bin/bash
# Holds the machine awake for as long as a given pid lives.
#
# The battery's own launcher takes sleep:idle, which covers an idle-triggered
# suspend. It does not cover the lid: logind handles that through a separate
# inhibitor type, handle-lid-switch, and this machine has HandleLidSwitch=suspend
# and HandleLidSwitchExternalPower=suspend. Closing the lid would therefore
# still have suspended a run. This takes every type that can stop one.
PID="$1"; WHY="$2"
exec systemd-inhibit --what=handle-lid-switch:sleep:idle:shutdown \
     --who="hold_awake for pid $PID" --why="$WHY" --mode=block \
     bash -c 'while kill -0 '"$PID"' 2>/dev/null; do sleep 30; done'
