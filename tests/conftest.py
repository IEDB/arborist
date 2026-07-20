import sys
from pathlib import Path

# Make the editable `protein_tree` package importable without a full install,
# so the offline regression tests run on any box (asahidake included).
SRC_PROTEIN = Path(__file__).parent.parent / 'src' / 'protein'
if str(SRC_PROTEIN) not in sys.path:
  sys.path.insert(0, str(SRC_PROTEIN))
