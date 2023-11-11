SELECT DISTINCT
  accession AS 'Accession',
  name AS 'Name',
  sequence AS 'Sequence',
  organism_id AS 'Organism ID'
FROM source
WHERE sequence IS NOT NULL
