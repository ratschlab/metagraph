#!/usr/bin/env bash
set -euo pipefail

# GitHub's Azure mirror can trickle downloads for tens of minutes without
# hitting APT's idle timeout: actions/runner-images#14823.
mirror_list=/etc/apt/apt-mirrors.txt
if [[ ! -f "$mirror_list" ]]; then
    exit 0
fi

hostname -f
sudo cp -p "$mirror_list" "${RUNNER_TEMP:?}/apt-mirrors.txt.original"
sudo tee "$mirror_list" > /dev/null <<'MIRRORS'
https://archive.ubuntu.com/ubuntu/	priority:1
https://security.ubuntu.com/ubuntu/	priority:2
MIRRORS
