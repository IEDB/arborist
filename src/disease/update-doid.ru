PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
PREFIX ONTIE: <https://ontology.iedb.org/ontology/ONTIE_>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

# move high-level nodes to 'additional disease by category'
DELETE { ?s rdfs:subClassOf obo:DOID_4 }
INSERT { ?s rdfs:subClassOf ONTIE:0003542 .
		 obo:DOID_417 rdfs:subClassOf obo:DOID_4 }
WHERE  {
	VALUES ?s { obo:DOID_7			# disease of anatomical entity
				obo:DOID_150		# disease of mental health
				obo:DOID_0014667	# disease of metabolism
				obo:DOID_630		# genetic disease
				obo:DOID_225		# syndrome
	}
} ;

# move asymptomatic terms under 'infection without disease'
INSERT { ?s rdfs:subClassOf ONTIE:0003424 }
WHERE  { ?s rdfs:label ?label .
		 FILTER STRSTARTS(?label, "asymptomatic")} ;

# move some to top-level diseases
DELETE { ?s rdfs:subClassOf ?parent }
INSERT { ?s rdfs:subClassOf obo:DOID_4 }
WHERE  {
	VALUES ?s { obo:DOID_1205 		# allergic disease
				obo:DOID_417 		# autoimmune disease
	}
	?s rdfs:subClassOf ?parent .
} ;

# change all IEDB alternative terms to true label
DELETE { ?s rdfs:label ?label }
INSERT { ?s rdfs:label ?alt ;
			oboInOwl:hasExactSynonym ?label }
WHERE  { ?s rdfs:label ?label ;
		    obo:OBI_9991118 ?alt . } ;


# move 'disease' under 'host status'
INSERT { obo:DOID_4 rdfs:subClassOf ONTIE:0003543 }
WHERE  {} ;


# Manual fix: autoimmune disease of skin and connective tissue should not be under musculoskeletal or rheumatic
DELETE { ?skinAndConn rdfs:subClassOf obo:DOID_0060032 ;
                      rdfs:subClassOf obo:DOID_1575 }
INSERT { ?skinAndConn rdfs:subClassOf obo:DOID_417 }
WHERE  { VALUES ?skinAndConn { obo:DOID_0060039 } } ;

# Manual fix: add 'reactive arthritis' under autoimmune
INSERT { ?reactArth rdfs:subClassOf obo:DOID_0060032 }
WHERE  { VALUES ?reactArth { obo:DOID_6196 } }
