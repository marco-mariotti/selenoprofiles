#!/usr/bin/env python
"""
"""
__author__  = "Francois Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


#import sys
#sys.path.insert(0, "/home/francisco/toolbox/utils")

import annotations.GO.Parsers.oboparser as obo
parserO = obo.Parser (open ("/home/francisco/toolbox/utils/annotations/GO/examples/data/gene_ontology_ext.obo"))
ontology = parserO.parse()

ontology.define_annotations_level()

# filter genes in associtation file that are "well" annotated
genes = {}
allowed_evidence = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'ISO', 'IC', 'TAS']
for annot in open("/home/francisco/toolbox/utils/annotations/GO/examples/data/gene_association.fb"):
    if annot.startswith ('!'): continue
    annot = annot.split('\t')
    if annot[6] not in allowed_evidence:
        continue
    genes.setdefault(annot [2], {})
    genes [annot [2]][annot[4]] = ontology.get_term_by_id (annot[4])

# propagate annotations
for g in genes:
    parents = []
    for annot in filter (lambda x: x.tags['level'] != 1, genes[g].values()):
        parents = list (set (parents + ontology.get_all_parents_ids (annot.id)))
    for dad in parents:
        genes[g][dad] = ontology.get_term_by_id (dad)

out = open ('/home/francisco/toolbox/utils/annotations/GO/examples/annot_fly.tab', 'w')
for gene in genes:
    for annot in genes [gene]:
        out.write (gene + '\t' + annot + '\n')

nam = 'GO:0070963'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0019051'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0070949'
ontology.draw_GO_plus_parents(nam)
nam = 'GO:0052040'
ontology.draw_GO_plus_parents(nam)




"""
1.  DB
Database from which annotated entry has been taken.
For the UniProtKB and Proteomes gene associaton files: UniProtKB
For the species-specific association files
(Arabidopsis, chicken, cow,
mouse, rat or zebrafish):
One of either: UniProtKB, ENSEMBL,
HINV (H-Invitational Database), TAIR, RefSeq or VEGA.
For the PDB association file:  PDB
Example: UniProtKB
2.  DB_Object_ID
A unique identifier in the database for the item being annotated.
Here: an accession number or identifier of the annotated protein
(or PDB entry for the gene_association.goa_pdb file)
For the UniProtKB, Human and Proteomes gene association files:
- a UniProtKB accession number or IPI identifier
For the IPI species-specific association files (Arabidopsis, chicken, cow,
mouse, rat or zebrafish):
- one of either UniProtKB, Ensembl, VEGA, HINV, TAIR or RefSeq peptide identifiers
For the PDB association file:
- a PDB identifier
Examples: O00165
ENSP00000241656
OTTDARP00000014036
HIT000018908
AT1G12760.2
NP_671756
117E
3.  DB_Object_Symbol
A (unique and valid) symbol (gene name) that corresponds to the DB_Object_ID.
An officially approved gene symbol will be used in this field when available.
Alternatively, other gene symbols or locus names are applied.
If no symbols are available, the identifier applied in column 2 will be used.
Examples: G6PC
CYB561
MGCQ309F3
C10H14ORF1
ENSBTAP00000000027
117E
4.  Qualifier
This column is used for flags that modify the interpretation of an
annotations.
If not null, then values in this field can equal: NOT, colocalizes_with, contributes_to,
NOT | contributes_to, NOT | colocalizes_with
Example: NOT
5.  GO ID
The GO identifier for the term attributed to the DB_Object_ID.
Example: GO:0005634
6.  DB:Reference
A single reference cited to support an annotations.
Where an annotations cannot reference a paper, this field will contain
a GO_REF identifier. See section 8 and
http://www.geneontology.org/doc/GO.references
for an explanation of the reference types used.
Examples: PMID:9058808
DOI:10.1046/j.1469-8137.2001.00150.x
GO_REF:0000002
GO_REF:0000020
GO_REF:0000004
GO_REF:0000003
GO_REF:0000019
GO_REF:0000023
GO_REF:0000024
7.  Evidence
One of either EXP, IMP, IC, IGI, IPI, ISS, IDA, IEP, IEA, TAS, NAS, NR, ND, ISO or RCA.
Example: TAS
  - Experimental Evidence Codes
    * EXP: Inferred from Experiment
    * IDA: Inferred from Direct Assay
    * IPI: Inferred from Physical Interaction
    * IMP: Inferred from Mutant Phenotype
    * IGI: Inferred from Genetic Interaction
    * IEP: Inferred from Expression Pattern
  - Computational Analysis Evidence Codes
    * ISS: Inferred from Sequence or Structural Similarity
    * ISO: Inferred from Sequence Orthology
    * ISA: Inferred from Sequence Alignment
    * ISM: Inferred from Sequence Model
    * IGC: Inferred from Genomic Context
    * RCA: inferred from Reviewed Computational Analysis
  - Author Statement Evidence Codes
    * TAS: Traceable Author Statement
    * NAS: Non-traceable Author Statement
  - Curator Statement Evidence Codes
    * IC: Inferred by Curator
    * ND: No biological Data available
  - Automatically-assigned Evidence Codes
    * IEA: Inferred from Electronic Annotation
  - Obsolete Evidence Codes
    * NR: Not Recorded
8.  With
An additional identifier to support annotationss using certain
evidence codes (including IEA, IPI, IGI, IC and ISS evidences).
Examples: UniProtKB:O00341
InterPro:IPROO1878
Ensembl:ENSG00000136141
GO:0000001
EC:3.1.22.1
9.  Aspect
One of the three ontologies, corresponding to the GO identifier applied.
P (biological process), F (molecular function) or C (cellular component).
Example: P
10. DB_Object_Name
Name of protein
The full UniProt protein name will be present here,
if available from UniProtKB. If a name cannot be added, this field
will be left empty.
Examples: Glucose-6-phosphatase
Cellular tumor antigen p53
Coatomer subunit beta
11. Synonym
Gene_symbol [or other text]
Alternative gene symbol(s), IPI identifier(s) and UniProtKB/Swiss-Prot identifiers are
provided pipe-separated, if available from UniProtKB. If none of these identifiers
have been supplied, the field will be left empty.
Example:  RNF20|BRE1A|IPI00690596|BRE1A_BOVIN
IPI00706050
MMP-16|IPI00689864
12. DB_Object_Type
What kind of entity is being annotated.
Here: protein (or protein_structure for the
gene_association.goa_pdb file).
Example: protein
13. Taxon_ID
Identifier for the species being annotated.
Example: taxon:9606
14. Date
The date of last annotations update in the format 'YYYYMMDD'
Example: 20050101
15. Assigned_By
Attribute describing the source of the annotations.  One of
either UniProtKB, AgBase, BHF-UCL, CGD, DictyBase, EcoCyc, EcoWiki, Ensembl, FlyBase, GDB, GeneDB_Spombe,
GeneDB_Pfal, GR (Gramene), HGNC, Human Protein Atlas, JCVI, IntAct, InterPro, LIFEdb, PAMGO_GAT, MGI, Reactome, RGD,
Roslin Institute, SGD, TAIR, TIGR, ZFIN, PINC (Proteome Inc.) or WormBase.
Example: UniProtKB
16. Annotation_Extension (N.B. Until annotations practices are finalised, this column will remain empty)
Contains cross references to other ontologies/databases that can be used to qualify or enhance the annotations.
The cross-reference is prefaced by an appropriate GO relationship; references to multiple ontologies
can be entered.
Example: part_of(CL:0000084)
occurs_in(GO:0009536)
has_input(CHEBI:15422)
has_output(CHEBI:16761)
has_participant(UniProtKB:Q08722)
part_of(CL:0000017)|part_of(MA:0000415)
17. Gene_Product_Form_ID
The unique identifier of a specific spliceform of the protein described in column 2 (DB_Object_ID)
Example:O43526-1
1.  DB
2.  DB_Object_ID
3.  DB_Object_Symbol
4.  Qualifier
5.  GO ID
6.  DB:Reference
7.  Evidence
8.  With
9.  Aspect
10. DB_Object_Name
11. Synonym
12. DB_Object_Type
13. Taxon_ID
14. Date
15. Assigned_By
16. Annotation_Extension
17. Gene_Product_Form_ID
genes.setdefault(annot [2], []).append (
        {'DB'                  : annot[0],
         'DB_Object_ID'        : annot[1],
         'DB_Object_Symbol'    : annot[2],
         'Qualifier'           : annot[3],
         'GO_id'               : annot[4],
         'DB:Reference'        : annot[5],
         'Evidence'            : annot[6],
         'With'                : annot[7],
         'Aspect'              : annot[8],
         'DB_Object_Name'      : annot[9],
         'Synonym'             : annot[10],
         'DB_Object_Type'      : annot[11],
         'Taxon_ID'            : annot[12],
         'Date'                : annot[13],
         'Assigned_By'         : annot[14],
         'Annotation_Extension': annot[15],
         'Gene_Product_Form_ID': annot[16].strip()}
        )
"""
