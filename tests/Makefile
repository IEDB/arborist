# Default: the fast, fully-offline suite. No bins, no MySQL, no network -- runs
# anywhere (asahidake included) in seconds. The bin/data-dependent e2e tests
# (test_assignments, test_leidos) self-skip when their prerequisites are absent.
test:
	@echo "Running offline test suite..."
	python3 -m pytest -q

# Full end-to-end on arborist-dev: generates a real run (needs blastp/mmseqs2 +
# ARC binaries in ../bin and build data), then asserts the golden outputs.
test-e2e:
	@echo "Running end-to-end suite (requires bins + build data)..."
	python3 ../src/protein/protein_tree/assign.py -b build -n 4
	python3 -m pytest -s test_assignments.py
