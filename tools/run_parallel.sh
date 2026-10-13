#!/usr/bin/env bash
# Run shell commands in parallel, wait for every process ID individually, and exit non-zero if
# any command failed (a bare `wait` would discard the individual exit statuses).
# Usage: tools/run_parallel.sh 'command 1' 'command 2' ...
set -u
pids=()
cmds=()
for cmd in "$@"; do
    bash -c "$cmd" &
    pids+=("$!")
    cmds+=("$cmd")
done
failed=0
for i in "${!pids[@]}"; do
    wait "${pids[$i]}"
    status=$?
    if [ "$status" -ne 0 ]; then
        echo "run_parallel: FAILED (exit $status): ${cmds[$i]}" >&2
        failed=$((failed + 1))
    fi
done
if [ "$failed" -ne 0 ]; then
    echo "run_parallel: $failed of ${#pids[@]} commands failed" >&2
    exit 1
fi
echo "run_parallel: all ${#pids[@]} commands succeeded"
