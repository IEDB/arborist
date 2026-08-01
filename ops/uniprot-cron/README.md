# UniProt release watch (arborist-dev)

Phase 1 of the UniProt-release → refresh → promote pipeline
(`/home/dmx2/reviews/uniprot-release-cron-feasibility.md` §10).

**This phase detects and notifies. It never runs `make`, never touches
pepmatch-curation, never writes a byte outside `/var/lib/arborist-uniprot-watch`.**

## What it does

Once a month, `watch-release.py` reads three independent signals:

| Signal | Why |
|---|---|
| `reldate.txt` (conditional GET on stored ETag) | 151 bytes; cheap 304 between releases |
| `rest.uniprot.org` `x-uniprot-release` header | what `select_proteome.py` actually queries |
| `reference_proteomes/RELEASE.metalink` `<version>` | says the *proteome set* landed, not just the KB |

`x-api-deployment-date` is ignored on purpose — the REST API redeploys without a
data release, so watching it produces false triggers.

All three must agree, and the release must have been out for ≥48 h
(`Last-Modified` on `reldate.txt`, or two sightings 48 h apart if that header is
missing). Only then does it write `RELEASE-CONFIRMED-<release>` into the state
dir and send a high-priority mail. That marker is the handoff point for the
refresh runner in a later phase; nothing consumes it yet.

Outcomes, all of which notify (silence must mean the timer is broken):

| Outcome | Exit | Effect |
|---|---|---|
| baseline seeded (first run) | 0 | adopts the current release, no marker |
| no new release | 0 | low-priority "timer ran" mail |
| sources disagree (mirror mid-sync) | 0 | heads-up mail, no marker |
| new release sighted, <48 h old | 0 | heads-up mail, no marker |
| new release confirmed | 0 | `RELEASE-CONFIRMED-<release>` + high-priority mail |
| poll failed | 70 | FAILED mail, `last_poll_error` in state, unit marked failed |
| not arborist-dev | 78 | refuses immediately |

## Install (run on arborist-dev, as root)

```sh
cd ~/lji/arborist                       # the arborist checkout on arborist-dev
git fetch origin && git checkout feat/uniprot-release-cron && git pull

sudo ops/uniprot-cron/install.sh "$PWD"
```

That installs `/opt/arborist-uniprot/bin/watch-release.py`, creates
`/var/lib/arborist-uniprot-watch/` and `/var/log/arborist-uniprot/`, writes both
systemd units with `ARBORIST_REPO` pointing at the checkout, and enables the
monthly timer. It refuses to run anywhere but arborist-dev.

Seed the baseline and prove the whole path works without waiting a month:

```sh
sudo /opt/arborist-uniprot/bin/watch-release.py --dry-run --no-notify   # reads only
sudo systemctl start arborist-uniprot-watch.service                     # real run + mail
systemctl status arborist-uniprot-watch.service
journalctl -u arborist-uniprot-watch.service -n 50 --no-pager
sudo cat /var/lib/arborist-uniprot-watch/state.json
systemctl list-timers arborist-uniprot-watch.timer
```

The first real run reports "baseline seeded at 2026_02" — that is correct, not a
missed release. The next distinct, settled release is the first one it announces.

## Operating it

```sh
# manual kick (same thing the timer does)
sudo systemctl start arborist-uniprot-watch.service

# change the schedule
sudo systemctl edit arborist-uniprot-watch.timer     # OnCalendar=*-*-05 06:30:00

# ship a code change after a git pull
sudo ops/uniprot-cron/install.sh "$PWD"

# stop it entirely
sudo systemctl disable --now arborist-uniprot-watch.timer
```

Optional ntfy push: uncomment `NTFY_URL` / `NTFY_TOPIC` in the service unit
(`NTFY_USER`/`NTFY_PASSWORD` if the server needs auth). Unset = email only.
Mail goes to `notify_email.DEFAULT_EMAIL_TO` (dmarrama@, rvita@); override with
`WATCH_EMAIL_TO` in the unit.

## Files

```
/opt/arborist-uniprot/bin/watch-release.py            the watcher
/var/lib/arborist-uniprot-watch/state.json            release memory (survives weekly_clean)
/var/lib/arborist-uniprot-watch/RELEASE-CONFIRMED-*   handoff marker for a later phase
/etc/systemd/system/arborist-uniprot-watch.{service,timer}
```

State lives under `/var/lib` deliberately: `make weekly_clean` wipes
`build/arborist/`, `cache/`, `current/` and would eat anything kept in the tree.

## Containment

Hostname gate in three places — `install.sh`, `watch-release.py`
(`exit 78` unless `hostname -s` is `arborist-dev`), and `ConditionHost=arborist-dev`
on both units. Never install this on asahidake or pepmatch-curation.

## Tests

`tests/test_watch_release.py`, offline, no network: `pytest tests/test_watch_release.py`.

## Not in this phase

`refresh-run.sh`, `diff-metrics.py`, `promote-to-prod.sh`, the build `flock`, k=3
index generation, and anything that writes to pepmatch-curation.
