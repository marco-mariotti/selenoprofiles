#!/usr/bin/python
"""
13 Nov 2010

Test building gene-ontology, annotations-extender and
displaying
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


import annotations.GO.Parsers.oboparser as obo
import annotations.GO.Parsers.annotparser as annotparser
from annotations.GO.annot import Annotation
from copy import copy


parserO = obo.Parser (open ("examples/data/gene_ontology_ext.obo"))
ontology = parserO.parse()

ontology.define_annotations_level()


parserA = annotparser.Parser (open ("examples/data/gene_association.fb"))
annotationss = parserA.parse(ontology)

genes = {}
allowed_evidence = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'ISO', 'IC', 'TAS']
for an in annotationss:
    if an.evidence_code._key not in allowed_evidence:
        continue
    genes.setdefault (an.db_object_id, []).append (an)

# propagate annotations
for g in genes:
    parents = []
    for go in map (lambda x: ontology.get_term_by_id (x.go_id), genes[g]):
        if go.tags['level']==1:
            continue
        parents = list (set (parents + ontology.get_all_parents_ids (go.id)))
    for dad in parents:
        item = copy (genes [g][0])
        item.assigned_by = 'propagation'
        item.go_id = dad
        genes[g].append (item)


nam = 'GO:0070963'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0019051'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0070949'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0052040'
ontology.draw_GO_plus_parents(nam)
