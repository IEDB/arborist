SELECT DISTINCT
  source_id AS 'Source ID',
  accession AS 'Accession',
  database AS 'Database',
  name AS 'Name',
  aliases AS 'Aliases',
  synonyms AS 'Synonyms',
  sequence AS 'Sequence',
  organism_id AS 'Organism ID'
FROM source
WHERE sequence IS NOT NULL
