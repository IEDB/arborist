SELECT DISTINCT
  epitope.epitope_id AS 'Epitope ID',
  object.mol1_seq AS 'Sequence',
  object.starting_position AS 'Starting Position',
  object.ending_position AS 'Ending Position',
  object.mol2_accession AS 'Source Accession',
  object.mol2_name AS 'Source Name',
  object.organism2_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.e_object_id = object.object_id
  AND object.mol1_seq IS NOT NULL
  AND object.object_sub_type IN 
      ('Peptide from protein', 'Discontinuous protein residues')

UNION

SELECT DISTINCT
  epitope.epitope_id AS 'Epitope ID',
  object.region AS 'Sequence',
  object.starting_position AS 'Starting Position',
  object.ending_position AS 'Ending Position',
  object.mol2_accession AS 'Source Accession',
  object.mol2_name AS 'Source Name',
  object.organism2_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.related_object_id = object.object_id
  AND object.region IS NOT NULL
  AND object.object_sub_type IN 
      ('Peptide from protein', 'Discontinuous protein residues')

UNION

SELECT DISTINCT
  epitope.epitope_id AS 'Epitope ID',
  object.region AS 'Sequence',
  object.starting_position AS 'Starting Position',
  object.ending_position AS 'Ending Position',
  object.mol2_accession AS 'Source Accession',
  object.mol2_name AS 'Source Name',
  object.organism2_id AS 'Organism ID'
FROM epitope, object
WHERE epitope.related_object_id = object.object_id
  AND object.mol1_seq IS NOT NULL
  AND object.object_sub_type IN ('Protein')
