#!/usr/bin/env bash
# Install the UniProt release watch + refresh runner on arborist-dev. Run as root.
#   sudo ops/uniprot-cron/install.sh [/path/to/arborist/repo]
# Idempotent: re-run after a git pull to ship a new watch-release.py.
set -euo pipefail

[[ "$(hostname -s)" == "arborist-dev" ]] || { echo "refusing: not arborist-dev (this is $(hostname -s))"; exit 78; }
[[ "$(id -u)" == "0" ]] || { echo "refusing: run as root (sudo)"; exit 77; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${1:-$(cd "$HERE/../.." && pwd)}"
[[ -f "$REPO/src/protein/protein_tree/notify_email.py" ]] || { echo "refusing: $REPO is not an arborist checkout"; exit 1; }

install -d -m 0755 /opt/arborist-uniprot/bin /var/lib/arborist-uniprot-watch /var/log/arborist-uniprot
install -m 0755 "$HERE/watch-release.py" /opt/arborist-uniprot/bin/watch-release.py
install -m 0755 "$HERE/refresh-run.sh" /opt/arborist-uniprot/bin/refresh-run.sh
install -m 0755 "$HERE/promote-to-prod.sh" /opt/arborist-uniprot/bin/promote-to-prod.sh
install -m 0644 "$HERE/README.md" /opt/arborist-uniprot/README.md

sed "s|^Environment=ARBORIST_REPO=.*|Environment=ARBORIST_REPO=$REPO|" \
  "$HERE/arborist-uniprot-watch.service" > /etc/systemd/system/arborist-uniprot-watch.service
install -m 0644 "$HERE/arborist-uniprot-watch.timer" /etc/systemd/system/arborist-uniprot-watch.timer
chmod 0644 /etc/systemd/system/arborist-uniprot-watch.service

systemctl daemon-reload
systemctl enable --now arborist-uniprot-watch.timer

echo "installed. repo=$REPO"
systemctl list-timers arborist-uniprot-watch.timer --no-pager
