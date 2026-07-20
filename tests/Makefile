
test:
	@echo "Testing Arborist process..."
	python3 ../src/protein/protein_tree/assign.py -b build -n 4
	python3 -m pytest -s test_assignments.py

# Fast, fully-offline regression for the taxon 1642 crash (bus brief #124).
# No MySQL, no organism rebuild, no network -- runs on any box in <1s.
test-1642:
	python3 -m pytest -q test_species_1642.py
