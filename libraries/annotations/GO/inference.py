#!/usr/bin/env python

"""Inference engine for drawing inferences on the Gene Ontology."""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["InferenceEngine", "Rules", "Query", "UNBOUND"]

from annotations.GO.ontology import GORelationship, GORelationshipFactory, Ontology
from collections import deque

UNBOUND = None

class Query(object):
    """A query on the Gene Ontology.
    
    Queries consist of three components: a subject term, a relation
    (e.g., ``is_a`` or ``part_of``) and an object term. Either one
    of the subject or the object terms can be unbound (represented
    by the `UNBOUND`_ constant, which maps to the ``None`` constant
    in Python).
    
    If both the subject and the object term is bound to a specific
    GO term, the result of the query will be either ``True`` or
    ``False``, depending on whether the relation can be inferred
    from the Gene Ontology using the standard inference rules or
    not.
    
    If the subject term is unbound, the result of the query will be
    a list of `GOTerm`_ instances that stand in the given relation
    with the given object term.
    
    If the object term is unbound, the result of the query will be
    a list of `GOTerm`_ instances that stand in the given relation
    with the given subject term.

    This class has three instance variables that you may modify:
    `subject_term`, `relation` and `object_term`. They refer to the
    subject term, the relation and the object term, respectively.
    """

    def __init__(self, ontology, subject_term=UNBOUND,
                 relation="is_a", object_term=UNBOUND):
        """Constructs a query on the Gene Ontology with the given
        subject and object terms and the given relation.
        
        :Parameters:
        - `ontology`: the ontology on which the query is formulated.
        - `subject_term`: the subject term as a `GOTerm` instance or its
          GO ID, or ``UNBOUND`` if the subject term is unbound in the
          query.
        - `relation`: the relation (a subclass of `GORelationship`).
        - `object_term`: the object term as a `GOTerm` instance or its
          GO ID, or ``UNBOUND`` if the object term is unbound in the
          query.
        """
        self._subject_term = None
        self._relation = None
        self._object_term = None

        self.ontology = ontology
        self.relation = relation
        self.subject_term = subject_term
        self.object_term = object_term

    def _get_relation(self):
        return self._relation
    def _set_relation(self, relation):
        if relation is UNBOUND:
            raise ValueError("relations may not be unbound")
        if not isinstance(relation, type) or not issubclass(relation, GORelationship):
            relation = GORelationshipFactory.from_name(relation)
        self._relation = relation
    relation = property(_get_relation, _set_relation,
                        doc="The relation part of the query")

    def _get_subject_term(self):
        return self._subject_term
    def _set_subject_term(self, subject_term):
        if subject_term is not UNBOUND:
            subject_term = self.ontology.ensure_term(subject_term)
        self._subject_term = subject_term
    subject_term = property(_get_subject_term, _set_subject_term,
                       doc="The subject term of the query")

    def _get_object_term(self):
        return self._object_term
    def _set_object_term(self, object_term):
        if object_term is not UNBOUND:
            object_term = self.ontology.ensure_term(object_term)
        self._object_term = object_term
    object_term = property(_get_object_term, _set_object_term,
                       doc="The object term of the query")


class InvalidQueryError(Exception):
    """Exception thrown when a query is invalid for some reason."""

    def __init__(self, query, reason):
        super.__init__(self, reason)
        self.query = query


def _ensure_relationship(rel):
    """Ensures that `rel` is a subclass of `GORelationship`. If it is a
    string, the name will be resolved using `GORelationshipFactory.from_name`."""
    if isinstance(rel, type) and issubclass(rel, GORelationship):
        return rel
    return GORelationshipFactory.from_name(rel)

class Rules(object):
    """Class that stores information about the rules of inference in the
    Gene Ontology.

    The summary of the `Gene Ontology inference rules`_ defines *rules*
    in the following format: "if we know that A is in a given relationship
    R1 with B, and B is in a given relationship R2 with C, then we can
    infer that A is in a given relationship R3 with C". Therefore, our
    rules will be stored in the form of triplets of R1, R2, and R3. For
    instance, the triplet ``("is_a", "is_a", "is_a")`` means the following
    rule: "if A ``is_a`` B and B ``is_a`` C, then A ``is_a`` C". Similarly,
    the triplet ``("regulates", "part_of", "regulates")`` means that
    "if A ``regulates`` B and B ``part_of`` C, then A ``regulates`` C".

    Besides rules, we also have *rule templates*. Rule templates differ
    from rules in the third element of the triplet, which may be either
    of the integer constants 1 or 2. 1 means that R3 is identical to R1,
    whatever R1 is, while 2 means that R3 is identical to R2. This enables
    us to encode all the five relationships where R1 is ``is_a`` as follows:
    ``("is_a", "any", 2)``, meaning that if A ``is_a`` B and B stands in
    an arbitrary relationship R2 with C, then A stands in the very same
    relationship with C.

    The rules and rule templates are stored in the `rules` member variable.
    There is also a class variable called `default_rules`, which contains
    the default inference rules of the Gene Ontology. If you don't supply
    a ruleset at construction time, the default rules will be used.

    Rules must be listen in the order they should be checked for; in other
    words, more specific rules (such as ``("positively_regulates", "is_a",
    "positively_regulates)``) should come before more general rules (such as
    ``("regulates", "inheritable", "regulates")``).
    """

    # TODO: check whether it's worth indexing the rules based on R1
    # TODO: is it worth refactoring stuff to a separate Rule class?

    default_rules = [
        ("regulates", "is_a", 1),
        ("regulates", "part_of", "regulates"),
        ("part_of", "part_of", "part_of"),
        ("part_of", "is_a", "part_of"),
        ("is_a", "any", 2),
    ]

    def __init__(self, rules=None):
        """Constructs a new ruleset with the given rules.

        If `rules` is ``None``, the default Gene Ontology inference rules
        will be used.
        """
        self.rules = []
        self._initialize(rules)

    def _initialize(self, rules=None):
        """Initializes the `rules` variable by transforming relationship
        names to `GORelationship` instances.

        If `rules` is ``None``, the default rules will be used.
        """
        if rules is None:
            rules = self.__class__.default_rules

        new_rules = []
        for rule in rules:
            target_rel = rule[2]
            if target_rel != 1 and target_rel != 2:
                target_rel = _ensure_relationship(target_rel)
            new_rule = (_ensure_relationship(rule[0]),
                        _ensure_relationship(rule[1]),
                        target_rel)
            new_rules.append(new_rule)
        self.rules = new_rules

    def apply(self, rel1, rel2):
        """Applies the inference rules to the two given relationships
        `rel1` and `rel2` and returns the inferred relationship or
        ``None`` if no inference can be drawn.
        """
        if rel1.object_term != rel2.subject_term:
            return None

        for rule in self.rules:
            if isinstance(rel1, rule[0]) and isinstance(rel2, rule[1]):
                if rule[2] == 1:
                    rel3 = rel1.__class__
                elif rule[2] == 2:
                    rel3 = rel2.__class__
                else:
                    rel3 = rule[2]
                return rel3(rel1.subject_term, rel2.object_term)
        return None

    def restrict_to_rules_relevant_for(self, relationship):
        """Constructs another `Rules`_ instance where only those inference
        rules will be kept that are relevant for inferring the given
        `relationship` (a subclass of `GORelationship`_).
        
        Returns the new `Rules`_ instance and the list of relationships
        that are relevant for inferring `relationship`.
        
        For instance, in order to infer `IsARelationship`_ instances, it is
        enough to consider other `IsARelationship`_ instances from the ontology.
        Therefore, the restricted ruleset will only contain rules involving
        `IsARelationship`_ relations, and the result will be a pair containing
        the new rules (as a `Rules` instance) and a singleton set containing
        `IsARelationship`_. In order to infer `PartOfRelationship`s,
        it is enough to consider `IsARelationship` and `PartOfRelationship`,
        therefore the restricted ruleset contains only the rules related
        to the inference of these two relations, plus a set containing these
        two relations.
        """
        if isinstance(relationship, GORelationship):
            relationship = relationship.__class__
        if not isinstance(relationship, type):
            raise TypeError("GORelationship or GORelationship instance expected")

        relevant_relationships = set([relationship])
        relevant_rule_idxs = set()
        # Find the transitive closure of relevant_relationships under the
        # current ruleset
        prev_count = 0
        while prev_count < len(relevant_relationships):
            prev_count = len(relevant_relationships)
            rels = relevant_relationships
            for idx, rule in enumerate(self.rules):
                # Have we seen this rule before? If we know that this rule
                # is relevant, we can skip it
                if idx in relevant_rule_idxs:
                    continue

                # Is it a rule of interest to us?
                # If we cannot use the rule to infer any relationship which
                # is not among our current relationships of interest or their
                # superclasses, continue
                if rule[2] == 1 or rule[2] == 2:
                    # This is a rule template
                    target_rel = rule[rule[2]-1]
                    if not any(issubclass(rel, target_rel) for rel in rels):
                        continue
                else:
                    if rule[2] not in rels:
                        continue

                # This is a new rule of interest.
                relevant_rule_idxs.add(idx)

                # We should be able to infer relationships standing on the left
                # hand side of the rule as well.
                if rule[2] == 1:
                    # Rule template. The first part of the rule doesn't matter
                    # as it is always equal to the conclusion, which is already
                    # in the set of relevant relationships.
                    relevant_relationships.add(rule[1])
                elif rule[2] == 2:
                    # Rule template. The second part of the rule doesn't matter
                    # as it is always equal to the conclusion, which is already
                    # in the set of relevant relationships.
                    relevant_relationships.add(rule[0])
                else:
                    # Exact rule, not a template. Both parts matter.
                    relevant_relationships.add(rule[0])
                    relevant_relationships.add(rule[1])

        return Rules(rule for idx, rule in enumerate(self.rules)
                     if idx in relevant_rule_idxs), relevant_relationships


class KnowledgeBase(object):
    """A collection of `GORelationship` instances.

    This class is really what it says on the tin: it is a collection of
    `GORelationship` instances. The only extra functionality offered by
    this class is that relationships are indexed by subject terms and
    object terms.
    """

    def __init__(self, relationships=None):
        """Constructs a knowledge base with the given relationships."""
        self._relationships = set()
        self._index_by_subject = {}
        self._index_by_object = {}

        if relationships:
            self.add_many(relationships)

    def __iter__(self):
        return iter(self.relationships)

    def __len__(self):
        return len(self.relationships)

    def add(self, relationship):
        """Adds a relationship to the knowledge base.
        
        Returns ``True`` if the relationship was not in the knowledge base, or
        ``False`` if the relationship was already known to the knowledge base.
        """
        num_rels = len(self._relationships)
        self._relationships.add(relationship)
        if len(self._relationships) == num_rels:
            return False

        if relationship.subject_term not in self._index_by_subject:
            self._index_by_subject[relationship.subject_term.id] = set()
        self._index_by_subject[relationship.subject_term.id].add(relationship)

        if relationship.object_term not in self._index_by_object:
            self._index_by_object[relationship.object_term.id] = set()
        self._index_by_object[relationship.object_term.id].add(relationship)

        return True

    def add_many(self, relationships):
        """Adds many relationships to the knowledge base.
        
        Returns those relationships which have not been added to the knowledge
        base before."""
        result = set(relationships)
        result.difference_update(self._relationships)
        if not result:
            return set()

        self._relationships.update(result)

        subjects = set(rel.subject_term.id for rel in result)
        objects = set(rel.object_term.id for rel in result)
        for term in subjects:
            if term not in self._index_by_subject:
                self._index_by_subject[term] = set()
        for term in objects:
            if term not in self._index_by_object:
                self._index_by_object[term] = set()

        for rel in result:
            self._index_by_subject[rel.subject_term.id].add(rel)
            self._index_by_object[rel.object_term.id].add(rel)

        return result

    def get_by_object_term(self, term):
        """Retrieves the set of relations where `term` is the object."""
        try:
            return set(self._index_by_object[term.id])
        except KeyError:
            return set()

    def get_by_subject_term(self, term):
        """Retrieves the set of relations where `term` is the subject."""
        try:
            return set(self._index_by_subject[term.id])
        except KeyError:
            return set()


class InferenceEngine(object):
    """Inference engine for the Gene Ontology."""

    def __init__(self, rules=None):
        """Constructs an inference engine that uses the given `rules` for
        inferring relations. `rules` must be an instance of `Rules` or
        ``None``, which means the default rules.
        """
        if rules is None:
            self.rules = Rules()
        else:
            self.rules = rules

    def solve(self, query, knowledge_base=None):
        """Solves the given query.
        
        Returns a generator that will iterate over relationships that match
        the query.

        :Parameters:
        - `query`: the query to be solved; an instance of `Query`.
        - `knowledge_base`: the knowledge base to be used for the inference. If
          ``None``, a new, empty knowledge base will be created. Specifying an
          existing knowledge base allows one to collect inferences regarding GO
          terms in a central place and re-use these inferences later; of
          course, such inferences are valid only if the ontology does not
          change in the meanwhile.
        """
        if knowledge_base is None:
            knowledge_base = KnowledgeBase()

        if query.subject_term is UNBOUND and query.object_term is UNBOUND:
            raise InvalidQueryError(query,
                    "subject and object are both unbound")

        if query.object_term is UNBOUND:
            return self.solve_unbound_object(query, knowledge_base)

        if query.subject_term is UNBOUND:
            raise NotImplementedError("unbound subject terms are not yet"
                                      "supported")

        return self.solve_bound(query, knowledge_base)

    def solve_unbound_object(self, query, knowledge_base=None):
        """Solves queries where the object is unbound and the subject is
        bound.

        Returns a generator that will iterate over relationships that match
        the query.

        You may call this method directly instead of `solve()`_ if you
        are sure that the query satisfies the preconditions of this method.
        """

        # TODO: cache the result of restrict_to_rules_relevant_for
        rules, rel_types = self.rules.restrict_to_rules_relevant_for(query.relation)
        rel_types = tuple(rel_types)

        subject_term = query.subject_term
        ontology = subject_term.ontology

        queued_term_ids = set([subject_term.id])
        terms_to_visit = deque([subject_term])

        if knowledge_base is None:
            knowledge_base = KnowledgeBase()

        while terms_to_visit:
            # Get an unvisited term
            term = terms_to_visit.popleft()

            # Get the relationships where this term is a subject
            rels = ontology.get_relationships(subject_term=term)
            if any(rel.subject_term == rel.object_term for rel in rels):
                raise ValueError

            # Filter to the relations of interest
            rels = [rel for rel in rels if isinstance(rel, rel_types)]

            # See if any of these new relations can be combined with the ones
            # in the knowledge base to infer new relations where the subject
            # is the one given originally
            if term is subject_term:
                new_rels = rels
            else:
                new_rels = []
                for rel2 in rels:
                    rels1 = knowledge_base.get_by_object_term(rel2.subject_term)
                    for rel1 in rels1:
                        rel3 = rules.apply(rel1, rel2)
                        if rel3 is not None:
                            new_rels.append(rel3)
            added_rels = knowledge_base.add_many(new_rels)

            # Iterate over the newly added relations and yield those which
            # we are interested in. Suppress those that we have seen before.
            for relation in added_rels:
                if isinstance(relation, query.relation):
                    yield relation

            # Add the newly discovered object terms to the list of nodes
            # we have to visit
            for rel in rels:
                if rel.object_term.id not in queued_term_ids:
                    terms_to_visit.append(rel.object_term)
                    queued_term_ids.add(rel.object_term.id)

    def solve_bound(self, query, knowledge_base=None):
        """Solves queries where the object and the subject are both bound.

        You may call this method directly instead of `solve()`_ if you
        are sure that the query satisfies the preconditions of this method.
        
        Queries of this type can be reformulated as follows: can it be inferred
        that the subject term is in the given relationship with the object
        term, given the inference rules and the ontology?
        
        Currently this is implemented by iterating over the results provided
        by `solve_unbound_object` and returning ``True`` as soon as we find
        the relationship in the conclusions drawn by `solve_unbound_object`.
        """
        expected = query.relation(query.subject_term, query.object_term)
        solutions = self.solve_unbound_object(query, knowledge_base)
        return any(solution == expected for solution in solutions)

