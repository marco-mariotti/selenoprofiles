"""
Some utility classes that are useful in `Bio.GO` but are unlikely
to be used in other parts of BioPython.

Maybe this module should be merged with `Bio.utils` eventually.
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"

__all__ = ["db_cursor_iterator", "pushback_iterator"]

def db_cursor_iterator(cursor, size=1):
    """Iterates over the results of the given database cursor.

    This enables us to use database cursors in ``for`` loops instead of
    a ``while`` loop, even for cursors that do not support iteration
    by design.
    """
    if hasattr(cursor, "__iter__"):
        for row in cursor:
            yield row
    while True:
        row = cursor.fetchone()
        if row is None:
            break
        yield row

# pylint:disable-msg=C0103
# C0103: invalid name
class pushback_iterator(object):
    """Constructs a pushback iterator from the given iterable.

    A pushback iterator is just like a regular iterable, but it
    has a `push_back` method that can be used to put some of
    the items back into the iterator temporarily.
    """

    def __init__(self, iterable):
        """Constructs a pushback iterator from the given iterable."""
        self.iterator = iter(iterable)
        self.buffer = []

    def __iter__(self):
        return self

    def next(self):
        """Returns the next item from the iterator, or the last
        one that was pushed back if there are items in the internal
        buffer."""
        try:
            return self.buffer.pop()
        except IndexError:
            return self.iterator.next()

    def push_back(self, item):
        """Pushes an item into the internal buffer so it will be
        returned when the user calls `next()` for the next time
        instead of the next element from the original iterable."""
        self.buffer.append(item)


def rewrite_db_query_markers(query, paramstyle):
    """Given an SQL query where the parameters are denoted by ``%s``
    markers, rewrites the query for the given parameter style.

    This method should be used for SQL queries when it is not known
    in advance which parameter style the database backend uses.

    It is assumed that ``%s`` occurs only unquoted (i.e. not within
    hard-coded string literals) in the query. It is also assumed that
    no other percent-style markers are present that could confuse
    Python's string substitution operator.

    :Parameters:
    - `paramstyle`: the desired parameter style. Must be one of:
      `format` (the query will not be touched at all),
      `qmark` (``%s`` will be replaced by ``?``) or
      `numeric` (``%s`` will be replaced by ``:1``, ``:2`` and so
      on, in this order).
    """

    if paramstyle == "format":
        return query

    if paramstyle == "qmark":
        return query.replace("%s", "?")

    if paramstyle == "numeric":
        n = query.count("%s")
        replacements = tuple(":%d" % i for i in xrange(1, n+1))
        return query % replacements

    raise NotImplementedError("unknown parameter style: %s" % paramstyle)
