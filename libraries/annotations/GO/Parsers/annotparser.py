#!/usr/bin/env python

"""
Parser for Gene Ontology annotations files

Usage example::

    import Bio.GO.Parsers.oboparser as obo
    import Bio.GO.Parsers.annotations as annotations

    parser = obo.Parser(open("gene_ontology.1_2.obo"))
    ontology = parser.parse()

    parser = annotations.Parser(open("gene_association.sgd"))
    annotationss = parser.parse(ontology)
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["ParseError", "Parser"]

from annotations.GO.annot import Annotation
from annotations.GO.utils import pushback_iterator

class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg, lineno = 1):
        Exception.__init__("%s near line %d" % (msg, lineno))
        self.lineno = lineno

class Parser(object):
    """The main attraction, the annotations file parser.
    
    If you want to create a parser that reads an annotations file, do this:

      >>> import Bio.GO.Parsers.annotations as annotations
      >>> parser = annotations.Parser(open("gene_associations.sgd"))

    Only the headers are read when creating the parser. You can
    access these right after construction as follows:

      >>> parser.headers["gaf-version"]
      ['1.0']

    To read the annotationss in the file, you must iterate over the
    parser as if it were a list. The iterator yields `Annotation`_
    objects. If you are not interested in the individual `Annotation`s
    and you only need an `Ontology` object, call the `parse()`
    method which will construct it for you.
    """

    # pylint:disable-msg=C0103
    def __init__(self, fp):
        """Creates an annotations parser that reads the given file-like object.
        """
        self.headers = {}

        if isinstance(fp, (str, unicode)):
            # not working..
            if fp.endswith(".gz"):
                # This is a gzipped file
                from gzip import GzipFile
                fp = GzipFile(fp)
            else:
                fp = open(fp)
        self._line_iterator = pushback_iterator(self._lines(fp))
        self._read_headers()

    @staticmethod
    def _lines(fp):
        """Iterates over the lines of a given file, removing
        comments and empty lines"""
        for line in fp:
            line = line.strip()
            if not line:
                continue
            yield line

    def _read_headers(self):
        """Reads the headers from the annotations file"""
        for line in self._line_iterator:
            if line[0] != '!':
                # We have reached the end of headers
                self._line_iterator.push_back(line)
                return

            parts = [part.rstrip() for part in line[1:].split(":", 1)]
            if len(parts) < 2:
                continue

            key, value = parts
            if key[0] == ' ':
                continue

            value = value.strip()
            try:
                self.headers[key].append(value)
            except KeyError:
                self.headers[key] = [value]

    # pylint:disable-msg=W0142
    # W0142: used * or ** magic
    def annotationss(self):
        """Iterates over the annotationss in this annotations file,
        yielding an `Annotation`_ object for each annotations."""
        for line in self._line_iterator:
            if line[0] == '!':
                continue
            parts = line.strip().split("\t")
            yield Annotation(*parts)

    def __iter__(self):
        return self.annotationss()

    def parse(self, ontology):
        """Parses the file handle given during construction time and
        returns an appropriately constructed `Annotations`_ instance.
        
        :Parameters:
            - `ontology`: the ontology being used to map term IDs to
              term names
        """
        for annotations in self:
            yield annotations
            # TODO
        pass
