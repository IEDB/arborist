

test:
	@echo "Testing Arborist process..."
	# python3 ../src/protein/protein_tree/select_proteome.py -b build/
	python3 ../src/protein/protein_tree/assign.py -b build -n 4 -t 108
	python3 -m pytest test_assignments.py
