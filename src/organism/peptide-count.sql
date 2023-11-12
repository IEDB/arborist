SELECT "Organism ID", COUNT('Sequence') AS 'Count'
FROM peptide
GROUP BY "Organism ID"
ORDER BY "Organism ID"
