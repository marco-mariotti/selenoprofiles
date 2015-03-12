#!/usr/bin/env python

"""Classes for handling annotationss of gene products to Gene Ontology
terms."""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["Annotation", "EvidenceCode"]

from annotations.Enum import Enum
#from annotations.GO.ontology import Aspect
from itertools import izip_longest

# pylint:disable-msg=W0232,R0903
# W0232: class has no __init__ method
# R0903: too few public methods
class EvidenceCode(Enum):
    """Evidence code of an annotations"""

    EXP = "Inferred from experiment"
    IDA = "Inferred from direct assay"
    IPI = "Inferred from physical interaction"
    IMP = "Inferred from mutant phenotype"
    IGI = "Inferred from genetic interaction"
    IEP = "Inferred from expression pattern"

    ISS = "Inferred from sequence or structural similarity"
    ISO = "Inferred from sequence orthology"
    ISA = "Inferred from sequence alignment"
    ISM = "Inferred from sequence model"
    IGC = "Inferred from genomic context"
    RCA = "Inferred from reviewed computational analysis"

    TAS = "Traceable author statement"
    NAS = "Non-traceable author statement"

    IC  = "Inferred by curator"
    ND  = "No biological data available"

    IEA = "Inferred from electronic annotations"

    NR  = "Not recorded"


class Annotation(object):
    """Class representing an entry from an annotations file.

    We support both GAF 1.0 and 2.0 annotations file formats. These files
    are simple flat files with 15 (GAF 1.0) or 17 (GAF 2.0) tab-delimited
    fields, some of which are required for all entries. Lines starting
    with an exclamation mark are silently ignored.

    An annotations consists of the following fields:

      - ``db``: refers to the database from which the identifier in
        ``db_object_id`` is drawn.

      - ``db_object_id``: a unique identifier from the database in
        ``db`` for the item being annotated.

      - ``db_object_symbol``: a (unique and valid) symbol to which
        ``db_object_id`` is matched. For instance, it can be the
        ORF name for an otherwise unnamed gene or protein. It should
        be a symbol that means something to a biologist whenever
        possible.

      - ``qualifiers``: flags that modify the interpretation of an
        annotations. Currently it can be one (or more) of ``NOT``,
        ``contributes_to`` and ``colocalizes_with``. As an annotations
        may have multiple qualifiers (or no qualifiers at all), this
        item is a tuple and it may be empty.

      - ``go_id``: the GO identifier for the term attributed to the
        ``db_object_id``.

      - ``db_references``: one or more unique identifiers for a
        *single* source cited as an authority for the attribution
        of the ``go_id`` to the ``db_object_id``. The general
        syntax is  ``DB:accession_number``. This should always
        be a tuple.

      - ``evidence_code``: the evidence code for the annotations.
        It is an instance of `EvidenceCode`_. Since it is mandatory
        in any GO annotations file, absent evidence codes will be
        replaced by `EvidenceCode.ND`.

      - ``froms``: required field for some evidence codes
        of this field. This item is a tuple of strings, each in the
        form of ``DB:accession_id`` or ``GO:GO_id``.

      - ``aspect``: refers to the namespace or ontology to which
        the given ``go_id`` belongs. It is an instance of
        `Aspect`_.

      - ``db_object_name``: name of gene or gene product. It may
        also be C{None}.

      - ``db_object_synonyms``: a tuple of alternative gene symbols
        that can be used in place of ``db_object_id``.

      - ``db_object_type``: a description of the type of gene product
        being annotated. It is usually one of ``protein_complex``,
        ``protein``, ``transcript``, ``ncRNA``, ``rRNA``, ``tRNA``,
        ``snRNA``, ``snoRNA``, ``gene_product`` or similar.

      - ``taxons``: a tuple of taxonomic identifiers.

      - ``date``: the date on which the annotations was made, in
        ``YYYYMMDD`` format.

      - ``assigned_by``: the database which made the annotations.

      - ``annotations_extensions``: a tuple of cross-references to
        other ontologies that can be used to qualify or enhance
        the annotations.

      - ``gene_product_form_id``: as the ``db_object_id`` entry
        *must* be a canonical entity (a gene or an abstract protein
        that has a 1:1 correspondence to a gene), this field allows
        the annotations of specific variants of that gene or gene
        product. It must be a standard 2-part global identifier
        such as ``UniProtKB:OK0206-2``.

    You may have noted that we used similar names to the ones used by
    the Gene Ontology Consortium in the GAF 2.0 Format Guide, but we
    tried to follow the Python style guide, and we also adopted a plural
    form for fields that may have cardinality greater than 1.
    """

    _fields = ("db", "db_object_id", "db_object_symbol", \
               "qualifiers", "go_id", "db_references", \
               "evidence_code", "froms", "aspect", \
               "db_object_name", "db_object_synonyms", \
               "db_object_type", "taxons", "date", \
               "assigned_by", "annotations_extensions", \
               "gene_product_form_id")
    __slots__ = _fields

    def __init__(self, *args, **kwds):
        """Constructs a new annotations object. Arguments should be given
        in the order they appear in a Gene Ontology annotations file.
        Keyword arguments may also be used and they take precedence over
        positional arguments.
        """
        for field, value in izip_longest(self._fields, args):
            setattr(self, field, value)
        for field, value in kwds.iteritems():
            setattr(self, field, value)
        self._tidy_fields()

    def _ensure_tuple(self, field_name):
        """Ensures that the field with the given name is a tuple"""
        value = getattr(self, field_name)
        if value is None or value == '':
            setattr(self, field_name, ())
        elif not isinstance(value, tuple):
            setattr(self, field_name, tuple(value.split("|")))

    # pylint:disable-msg=W0201,E0203,E1101
    # W0201: attribute defined outside __init__
    # E0203: access to member before its definition
    # E1101: class has no 'from_value' member
    def _tidy_fields(self):
        """Tidies up the fields in this annotations by casting them to the
        required format."""
        self._ensure_tuple("qualifiers")
        self._ensure_tuple("db_references")
        self._ensure_tuple("db_object_synonyms")
        self._ensure_tuple("froms")
        self._ensure_tuple("taxons")
        self._ensure_tuple("annotations_extensions")

        if self.evidence_code is None:
            self.evidence_code = EvidenceCode.ND
        elif not isinstance(self.evidence_code, EvidenceCode):
            self.evidence_code = EvidenceCode.from_name(self.evidence_code)

        #if self.aspect is not None and not isinstance(self.aspect, Aspect):
        #    self.aspect = Aspect.from_name(self.aspect)

    def __repr__(self):
        """Returns a Python representation of this object"""
        fields = ["%s=%r" % (field, getattr(self, field)) \
                  for field in self._fields]
        return "%s(%s)" % (self.__class__.__name__, ", ".join(fields))



