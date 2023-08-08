-- Get a list of all NCBI taxa used in IEDB data
SELECT ncbi_organism_id AS taxon_id
FROM ncbi_organism_include
WHERE ncbi_organism_id < 10000000
ORDER BY taxon_id
