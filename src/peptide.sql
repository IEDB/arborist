SELECT DISTINCT
  object.mol1_seq AS 'Sequence',
  object.mol2_name AS 'Source Name',
  object.mol2_accession AS 'Accession',
  object.organism2_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.e_object_id = object.object_id
  AND object.mol1_seq IS NOT NULL
  AND object.object_sub_type IN 
      ('Peptide from protein', 'Discontinuous protein residues')

UNION

SELECT DISTINCT
  object.region AS 'Sequence',
  object.mol2_name AS 'Source Name',
  object.mol2_accession AS 'Accession',
  object.organism2_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.related_object_id = object.object_id
  AND object.region IS NOT NULL
  AND object.object_sub_type IN 
      ('Peptide from protein', 'Discontinuous protein residues')

UNION

SELECT DISTINCT
  object.mol1_seq AS 'Sequence',
  object.mol1_name AS 'Source Name',
  object.mol1_accession AS 'Accession',
  object.organism_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.related_object_id = object.object_id
  AND object.mol1_seq IS NOT NULL
  AND object.object_sub_type IN ('Protein')
