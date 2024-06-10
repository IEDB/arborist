SELECT DISTINCT
  epitope.epitope_id AS 'Epitope ID',
  COALESCE(object.mol1_seq, object.region) AS 'Sequence',
  object.starting_position AS 'Starting Position',
  object.ending_position AS 'Ending Position',
  object.mol2_accession AS 'Source Accession',
  object.organism2_id AS 'Organism ID',
  object.organism2_name AS 'Organism Name'
FROM
  object
JOIN
  epitope ON epitope.e_object_id = object.object_id
WHERE
  object.object_sub_type IN ('Peptide from protein', 'Discontinuous protein residues')

UNION

SELECT DISTINCT 
  epitope.epitope_id AS 'Epitope ID',
  COALESCE(object.mol1_seq, object.region) AS 'Sequence',
  object.starting_position AS 'Starting Position',
  object.ending_position AS 'Ending Position',
  object.mol2_accession AS 'Source Accession',
  object.organism2_id AS 'Organism ID',
  object.organism2_name AS 'Organism Name'
FROM
  object
JOIN
  epitope ON epitope.related_object_id = object.object_id
WHERE
  object.object_sub_type IN ('Peptide from protein', 'Discontinuous protein residues')