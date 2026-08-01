# UniProt release watch (arborist-dev)

The UniProt-release → refresh → promote pipeline
(`/home/dmx2/reviews/uniprot-release-cron-feasibility.md` §10), as far as it is built.

**Nothing here promotes anything.** `watch-release.py` detects releases and mails
you; `refresh-run.sh` rebuilds arborist-dev. Neither writes a byte to
pepmatch-curation.

## What it does

Once a week, `watch-release.py` reads three independent signals:

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
cd /mnt/data/arborist
git checkout main && git pull

sudo ops/uniprot-cron/install.sh "$PWD"
```

That installs `/opt/arborist-uniprot/bin/watch-release.py`, creates
`/var/lib/arborist-uniprot-watch/` and `/var/log/arborist-uniprot/`, writes both
systemd units with `ARBORIST_REPO` pointing at the checkout, and enables the
weekly timer. It refuses to run anywhere but arborist-dev.

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
(`exit 78` unless `hostname -s` is `arborist-dev`), and `ConditionHost=arborist-dev*`
on both units (glob, because the box reports the FQDN `arborist-dev.lji.org`).

## The refresh (`refresh-run.sh`)

Run by hand after the watch reports a confirmed release. Nothing schedules it.

```sh
sudo /opt/arborist-uniprot/bin/refresh-run.sh            # DATESTAMP auto-resolved
sudo /opt/arborist-uniprot/bin/refresh-run.sh 20260726   # or pin the snapshot
```

It takes `flock -n /var/lock/arborist-build.lock` so a weekly build can never
overlap, wipes `/mnt/data/arborist/build/species/` (literal path, no variable),
then runs:

```
make weekly_clean
make deps iedb ncbitaxon organism proteome protein
```

`proteome` is the point — `make weekly` omits it, so `arborist.sh` alone would
reassign against stale FASTA and never call `select_proteome.py`. The run goes
all the way through `protein` so arborist-dev is a usable snapshot and the
assignment numbers exist for a later gate. It promotes nothing.

`DATESTAMP` is the most recent Sunday, then walks back up to 4 Sundays until an
`iedb_query_<stamp>` database exists on `web02.internal.iedb.org`, logging a
warning if it had to fall back. It aborts if `/mnt/data` is under 200 G free,
because `assign.py` reuses any per-species output with size > 0 and a
disk-death mid-run leaves truncated files the next run silently trusts.

Afterwards, `build/reports/release-stamp.tsv` records which UniProt release,
IEDB snapshot and git SHA produced the tree — `.pepidx` carries no release
field, so that file is the only answer to "what is this?".

Expect hours: the species wipe means every proteome is re-fetched from UniProt.

## Tests

Offline, no network, no arborist-dev:

```sh
pytest tests/test_watch_release.py tests/test_refresh_run.py
```

## Not built yet

`diff-metrics.py` (the numbers gate), `promote-to-prod.sh`, `flock` on the
weekly's own schedule, k=3 index generation, and anything that writes to
pepmatch-curation.
