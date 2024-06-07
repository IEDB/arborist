SELECT DISTINCT
  accession AS 'Source Accession',
  database AS 'Database',
  name AS 'Name',
  aliases AS 'Aliases',
  synonyms AS 'Synonyms',
  sequence AS 'Sequence',
  organism_id AS 'Organism ID',
  organism_name AS 'Organism Name',
  iri as 'IRI'
FROM source
