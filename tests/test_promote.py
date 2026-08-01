"""Guards on the only script that writes outside arborist-dev."""

import re
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent.parent / 'ops' / 'uniprot-cron' / 'promote-to-prod.sh'
SOURCE = SCRIPT.read_text()


def run(*args):
  return subprocess.run(['bash', str(SCRIPT), *args], capture_output=True, text=True)


def test_refuses_without_an_explicit_mode():
  """No accidental promotion from a bare invocation."""
  result = run()
  assert result.returncode != 0
  assert '--dry-run' in result.stderr and '--confirm' in result.stderr


def test_refuses_off_the_build_host():
  result = run('--confirm')
  assert result.returncode != 0
  assert 'arborist-dev.lji.org' in result.stderr


def test_host_gate_precedes_every_remote_write():
  """The gate must come before ssh, rsync and the snapshot, in source order."""
  gate = SOURCE.index('BUILD_HOST_FQDN=')
  for command in ('rsync ', 'ssh "${SSH_OPTS', 'cp -al'):
    assert gate < SOURCE.index(command), f'{command} appears before the host gate'


def test_gate_requires_clean_diff_or_explicit_approval():
  assert 'last_gate_exit' in SOURCE
  assert 'APPROVED-$release' in SOURCE
  # the refusal has to tell you how to approve, or the token is folklore
  assert 'mv $STATE_DIR/HOLD-$release' in SOURCE


def test_prod_is_a_full_mirror_not_a_delta():
  """Dev and prod must match exactly; syncing only changed species would drift."""
  prod_sync = next(line for line in SOURCE.splitlines()
                   if 'PROD_HOST:$PROD_SPECIES/' in line and line.strip().startswith('rsync'))
  assert '--delete' in prod_sync
  assert '$SPECIES_DIR/' in prod_sync


def test_snapshot_is_hardlinked_and_dropped_on_success():
  """A copy would need 83G of the 104G free on prod; hardlinks need none."""
  assert 'cp -al' in SOURCE
  assert 'cp -a ' not in SOURCE
  assert re.search(r'rm -rf .\$snapshot', SOURCE)
  # dropped only after both verifications
  assert SOURCE.index('verifying pepmatch') < SOURCE.index('dropping snapshot')


def test_pepmatch_sync_cannot_leave_a_torn_index():
  """A reader mmapping a half-written .pepidx is the failure mode to avoid."""
  pepmatch_block = SOURCE[SOURCE.index('pepmatch_rsync=('):SOURCE.index('# --- verify')]
  assert '--temp-dir=' in pepmatch_block
  assert '--inplace' not in SOURCE
  assert '--partial ' not in SOURCE  # --partial-dir on the prod tree is fine


def test_pepmatch_uses_the_existing_key():
  """Reuse dmarrama's key rather than minting a second credential."""
  assert 'PEPMATCH_KEY=/home/dmarrama/.ssh/pepmatch_curation_rsa' in SOURCE
  assert 'IdentitiesOnly=yes' in SOURCE
  assert 'BatchMode=yes' in SOURCE


def test_indexes_are_flattened_for_pepmatch():
  """Destination is one flat directory: <taxon>_<k>mers.pepidx."""
  assert '${taxon}_${kmer}.pepidx' in SOURCE
  assert 'proteome_*mers.pepidx' in SOURCE


def test_failure_after_prod_sync_keeps_the_rollback_point():
  for message in ('prod rsync failed', 'pepmatch rsync failed', 'prod differs after sync'):
    line = next(l for l in SOURCE.splitlines() if message in l)
    assert 'snapshot' in line, f'{message} does not mention the rollback point'


@pytest.mark.parametrize('taxon', [9606, 10090, 10116, 9612, 9796, 9823, 9913, 9986])
def test_three_mer_species_list(taxon):
  from protein_tree.preprocess_3mers import THREE_MER_SPECIES
  assert taxon in THREE_MER_SPECIES


def test_three_mer_list_is_exactly_the_eight_served():
  from protein_tree.preprocess_3mers import THREE_MER_SPECIES, K
  assert len(THREE_MER_SPECIES) == 8
  assert K == 3


def test_three_mer_build_skips_current_indexes(tmp_path):
  from protein_tree.preprocess_3mers import needs_build, index_path
  species = tmp_path / '9606'
  species.mkdir()
  fasta = species / 'proteome.fasta'
  fasta.write_text('>a\nMK\n')

  assert needs_build(species) is True

  index = index_path(species)
  index.write_text('x')
  import os
  os.utime(index, (fasta.stat().st_mtime + 10, fasta.stat().st_mtime + 10))
  assert needs_build(species) is False
  assert needs_build(species, force=True) is True


def test_three_mer_build_skips_empty_proteomes(tmp_path):
  """An EMPTY species has nothing to index; that is not an error."""
  from protein_tree.preprocess_3mers import needs_build
  species = tmp_path / '1234'
  species.mkdir()
  (species / 'proteome.fasta').write_text('')
  assert needs_build(species) is False
