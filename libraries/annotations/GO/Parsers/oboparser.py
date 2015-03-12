#!/usr/bin/env python
"""
A very simple and not 100% compliant parser for the OBO file format

This parser is supplied "as is". It is not an official parser, it
might refuse to parse perfectly valid OBO files, or it might parse
perfectly invalid OBO files; on the other hand, it can parse the
official Gene Ontology OBO files and the ones created by OBO-Edit,
so it isn't too bad :)

Usage example::

    import Bio.GO.Parsers.oboparser as obo
    parser = obo.Parser(open("gene_ontology.1_2.obo"))
    ontology = parser.parse()
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"

__all__ = ["ParseError", "Stanza", "Parser", "Value"]


from cStringIO import StringIO
from warnings import warn

import re
import tokenize

from annotations.GO.ontology import GeneOntologyNX, GOTerm, Ontology, \
     IsARelationship, PartOfRelationship, \
     RegulatesRelationship, \
     PositivelyRegulatesRelationship, \
     NegativelyRegulatesRelationship

class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg, lineno = 1):
        Exception.__init__("%s near line %d" % (msg, lineno))
        self.lineno = lineno


class Value(object):
    """Class representing a value and its modifiers in the OBO file

    This class has two member variables. `value` is the value itself,
    `modifiers` are the corresponding modifiers in a tuple. Currently
    the modifiers are not parsed in any way, but this might change in
    the future.
    """

    __slots__ = ["value", "modifiers"]

    def __init__(self, value, modifiers=()):
        """Creates a new value"""
        self.value = str(value)
        if modifiers:
            self.modifiers = tuple(modifiers)
        else:
            self.modifiers = None

    def __eq__(self, other):
        """Tests whether this `Value` instance is equal to another one.
        
        Two `Value` instances are equal if they have equal value with the
        same modifiers in the same order.
        """
        return self.value == other.value and self.modifiers == other.modifiers

    def __str__(self):
        """Returns the value itself (without modifiers)"""
        return str(self.value)

    def __repr__(self):
        """Returns a Python representation of this object"""
        return "%s(%r, %r)" % (self.__class__.__name__, \
                self.value, self.modifiers)


class Stanza(object):
    """Class representing an OBO stanza.

    An OBO stanza looks like this::

      [name]
      tag: value
      tag: value
      tag: value

    Values may optionally have modifiers, see the OBO specification
    for more details. This class stores the stanza name in the
    `name` member variable and the tags and values in a Python
    dict called `tags`. Given a valid stanza, you can do stuff like
    this:

      >>> stanza.name
      "Term"
      >>> print stanza.tags["id"]
      ['GO:0015036']
      >>> print stanza.tags["name"]
      ['disulfide oxidoreductase activity']

    Note that the `tags` dict contains lists associated to each
    tag name. This is because theoretically there could be more than
    a single value associated to a tag in the OBO file format.
    """

    __slots__ = ["name", "tags"]

    def __init__(self, name, tags=None):
        """Creates a new stanza with the given name and the given
        tags (which must be a dict)"""
        self.name = name
        if tags:
            self.tags = dict(tags)
        else:
            self.tags = dict()

    def __repr__(self):
        """Returns a Python representation of this object"""
        return "%s(%r, %r)" % (self.__class__.__name__, \
                self.name, self.tags)

    def add_tag(self, name, value):
        """Adds a tag-value pair to this stanza. If the tag name already
        exists, the value will be appended to the value list of that tag.
        """
        try:
            self.tags[name].append(value)
        except KeyError:
            self.tags[name] = [value]

    def to_term(self, term_factory=GOTerm):
        """Converts this stanza to an instance of the given term class.

        The stanza must contain a term; in other words, the stanza name must
        be equal to ``"Term"``. This is not checked in this method.

        :Parameters:
        - `term_factory`: the term factory to be used. It is safe to leave
          it at its default value unless you want to use a custom term class.
          The signature of this factory function should be identical to that
          of the constructor of `GOTerm`.
        """
        identifier = str(self.tags["id"][0])
        name = str(self.tags.get("name", [identifier])[0])
        aliases = [str(go_id) for go_id in self.tags.get("alt_id", [])]
        return term_factory(identifier, name, aliases, self.tags)


class Parser(object):
    """The main attraction, the OBO parser.
    
    If you want to create a parser that reads an OBO file, do this:

      >>> import annotations.GO.Parsers.oboparser as obo
      >>> parser = obo.Parser(open("gene_ontology.1_2.obo"))

    Only the headers are read when creating the parser. You can
    access these right after construction as follows:

      >>> parser.headers["format-version"]
      ['1.2']

    To read the stanzas in the file, you must iterate over the
    parser as if it were a list. The iterator yields `Stanza`
    objects. If you are not interested in the individual `Stanza`s
    and you only need an `Ontology` object, call the `parse()`
    method which will construct it for you.
    """

    def __init__(self, fp):
        """Creates an OBO parser that reads the given file-like object.
        """
        self._line_buffer = []

        if isinstance(fp, (str, unicode)):
            fp = open(fp)
        self.fp = fp
        self.tag_value_pair_re = re.compile(r"\s*(?P<tag>[^:]+):\s*(?P<value>.*)")
        self.stanza_name_re = re.compile(r"\[(?P<name>[^]]*)\]")
        self.lineno = 0
        self._read_headers()

    def _lines(self):
        """Iterates over the lines of the file, removing
        comments and trailing newlines and merging multi-line
        tag-value pairs into a single line"""
        while self._line_buffer:
            yield self._line_buffer.pop(-1)

        while True:
            self.lineno += 1
            line = self.fp.readline()
            if not line:
                break

            line = line.strip()
            if not line:
                # This is an empty line
                yield line
                continue

            if line[0] == '!':
                # This is a comment line, so it can be ignored
                continue

            if line[-1] == '\\':
                # This line is continued in the next line
                lines = [line[:-1]]
                finished = False
                while not finished:
                    self.lineno += 1
                    line = self.fp.readline()
                    if not line:
                        # End of file. Well, this should not happen in a valid
                        # OBO file, but we are not that strict
                        finished = True
                        break

                    if line[0] == '!':
                        # A comment line in the middle of a continuation
                        continue

                    line = line.strip()
                    if line and line[-1] == '\\':
                        # Continuation follows in the next line
                        lines.append(line[:-1])
                    else:
                        # This is the final line of this block
                        lines.append(line)
                        finished = True
                line = " ".join(lines)
            else:
                # No line continuation
                try:
                    # Search for a trailing comment
                    # The OBO specification is a bit vague here as it does not
                    # specify whether trailing comments can contain ! or not.
                    # We assume that they cannot - this works for the Gene
                    # Ontology for the time being.
                    comment_char = line.rindex("!")
                    line = line[0:comment_char].strip()
                except ValueError:
                    # No comment, fine
                    pass

            yield line

    def _parse_tag_value_pair(self, line):
        """Parses a single line consisting of a tag-value pair
        and optional modifiers. Returns the tag name and the
        value as a `Value` object, or ``False`` if the line does
        not contain a tag-value pair."""
        match = self.tag_value_pair_re.match(line)
        if not match:
            return False

        tag, value_and_mod = match.group("tag"), match.group("value")

        # If the value starts with a quotation mark, we parse it as a
        # Python string -- luckily this is the same as an OBO string
        if value_and_mod and value_and_mod[0] == '"':
            g = tokenize.generate_tokens(StringIO(value_and_mod).readline)
            for toknum, tokval, _, (erow, ecol), _ in g:
                if toknum == tokenize.STRING:
                    value = eval(tokval)
                    mod = (value_and_mod[ecol:].strip(), )
                    break
                raise ParseError("cannot parse string literal", self.lineno)
            return tag, Value(value, mod)

        return tag, Value(value_and_mod, None)

    def _read_headers(self):
        """Reads the headers from the OBO file"""
        self.headers = {}
        for line in self._lines():
            if not line or line[0] == '[':
                # We have reached the end of headers
                # Push back the current line to the buffer as we
                # will need it later once again
                self._line_buffer.append(line)
                return
            key, value = self._parse_tag_value_pair(line)
            try:
                self.headers[key].append(value.value)
            except KeyError:
                self.headers[key] = [value.value]

    def stanzas(self):
        """Iterates over the stanzas in this OBO file,
        yielding a `Stanza` object for each stanza."""
        stanza = None
        stanza_name_re = self.stanza_name_re

        for line in self._lines():
            if not line:
                continue

            # Do we have a stanza name in this line?
            match = stanza_name_re.match(line)
            if match:
                # Yes, yield the current stanza and start a new one
                if stanza:
                    yield stanza
                stanza = Stanza(match.group("name"))
            else:
                # No, the line contains a tag-value pair
                tag, value = self._parse_tag_value_pair(line)
                stanza.add_tag(tag, value)

        # Yield the last stanza (if any)
        if stanza:
            yield stanza

    def __iter__(self):
        return self.stanzas()

    def parse(self, load_obsolete=False, ontology_factory=GeneOntologyNX, \
              term_factory=Stanza.to_term):
        """Parses the file handle given during construction time and
        returns an appropriately constructed `Ontology` instance.
        
        :Parameters:
            - `load_obsolete`: whether to load obsolete entries from the
              ontology file.
            - `ontology_factory`: a factory that generates an empty instance
              of `Ontology`. The default is safe, you have to override
              its value only if you want to use a custom ontology class.
            - `term_factory`: a factory that generates an appropriate `Term`
              instance from a `Stanza` instance. The default is safe, you have
              to override its value only if you want to use a custom term
              class in the ontology.
        """
        ontology = ontology_factory()
        relationships = []

        if load_obsolete:
            stanza_iter = (stanza for stanza in self if stanza.name == "Term")
        else:
            true_value = Value("true")
            stanza_iter = (stanza for stanza in self \
                           if stanza.name == "Term" and \
                           true_value not in stanza.tags.get("is_obsolete", []))

        # Dict to map relationship types to GORelationship classes
        rel_types = {
            "is_a": IsARelationship,
            "part_of": PartOfRelationship,
            "regulates": RegulatesRelationship,
            "positively_regulates": PositivelyRegulatesRelationship,
            "negatively_regulates": NegativelyRegulatesRelationship
        }

        # Add the terms
        for stanza in stanza_iter:
            term = term_factory(stanza)
            ontology.add_term(term)

            rel = IsARelationship
            for ancestor in stanza.tags.get("is_a", []):
                relationships.append((term, rel, ancestor.value))

            for value in stanza.tags.get("relationship", []):
                rel_type, value = value.value.split(None, 1)
                rel = (term, rel_types.get(rel_type, None), value)
                if rel[1]:
                    relationships.append(rel)
                ## my comment, too verbose
                #else:
                #    warn("ignoring unknown relationship: %r %s %r" % rel)

        # Add the relationships
        for subject_term, rel, ancestor in relationships:
            # We have to look up the ancestor in the ontology
            object_term = ontology.get_term_by_id(ancestor)
            ontology.add_relationship(subject_term=subject_term, \
                                      object_term=object_term, \
                                      relationship=rel)

        return ontology
