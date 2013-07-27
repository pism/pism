# Copyright (c) 2007, 2008, 2009, 2010, 2011, 2012  Andrey Golovizin
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""Miscellaneous small utils."""


from functools import wraps
from collections import Sequence, MutableMapping
from types import GeneratorType


def memoize(f):
    memory = {}
    @wraps(f)
    def new_f(*args):
        if args not in memory:
            memory[args] = f(*args)
        return memory[args]
    return new_f


class CaseInsensitiveDict(MutableMapping):
    """A dict with case-insensitive lookup.

    >>> d = CaseInsensitiveDict(TesT='passed')
    >>> d
    CaseInsensitiveDict({'TesT': 'passed'})
    >>> print d['TesT']
    passed
    >>> print d['test']
    passed
    >>> print d['Test']
    passed
    >>> print d.get('test')
    passed
    >>> print d.get('Test')
    passed

    >>> d['Test'] = 'passed again'
    >>> print d['test']
    passed again
    >>> 'test' in d
    True
    >>> print d.keys()
    ['Test']
    >>> for key in d:
    ...     print key
    Test
    >>> for key, value in d.iteritems():
    ...     print key, value
    Test passed again
    >>> bool(d)
    True
    >>> len(d)
    1

    >>> del d['test']
    >>> len(d)
    0
    >>> bool(d)
    False
    >>> 'test' in d
    False
    >>> 'Test' in d
    False
    >>> print d.get('test')
    None
    >>> print d.get('Test')
    None
    >>> print d.get('Test', 'failed')
    failed

    >>> CaseInsensitiveDict(
    ...     (key, value) for key, value in [('a', 'b')]
    ... )
    CaseInsensitiveDict({'a': 'b'})

    """

    def __init__(self, *args, **kwargs):
        self._dict = dict(*args, **kwargs)
        self._keys = dict((key.lower(), key) for key in self)

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._dict)

    def __setitem__(self, key, value):
        """To implement lowercase keys."""
        key_lower = key.lower()
        try:
            existing_key = self._keys[key_lower]
        except KeyError:
            pass
        else:
            del self._dict[existing_key]
        self._dict[key] = value
        self._keys[key_lower] = key

    def __getitem__(self, key):
        existing_key = self._keys[key.lower()]
        return self._dict[existing_key]

    def __delitem__(self, key):
        key_lower = key.lower()
        existing_key = self._keys[key_lower]
        del self._dict[existing_key]
        del self._keys[key_lower]

    def __deepcopy__(self, memo):
        from copy import deepcopy
        return CaseInsensitiveDict(
            (key, deepcopy(value, memo)) for key, value in self.iteritems()
        )

    def __contains__(self, key):
        return key.lower() in self._keys

    def __repr__(self):
        """A caselessDict version of __repr__ """
        return '{0}({1})'.format(
            type(self).__name__, self._dict.__repr__()
        )


class CaseInsensitiveDefaultDict(CaseInsensitiveDict):
    """CaseInseisitiveDict with default factory, like collections.defaultdict
    
    >>> d = CaseInsensitiveDefaultDict(int)
    >>> d['a']
    0
    >>> d['a'] += 1
    >>> d['a']
    1
    >>> d['A']
    1
    >>> d['a'] = 3
    >>> d['a']
    3
    >>> d['B'] += 10
    >>> d['b']
    10

    """
    def __init__(self, default_factory):
        super(CaseInsensitiveDefaultDict, self).__init__()
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return super(CaseInsensitiveDefaultDict, self).__getitem__(key)
        except KeyError:
            return self.default_factory()


class OrderedCaseInsensitiveDict(CaseInsensitiveDict):
    """ An (incomplete) ordered case-insensitive dict.

    >>> d = OrderedCaseInsensitiveDict([
    ...     ('uno', 1),
    ...     ('dos', 2),
    ...     ('tres', 3),
    ... ])
    >>> d.keys()
    ['uno', 'dos', 'tres']
    >>> d.items()
    [('uno', 1), ('dos', 2), ('tres', 3)]
    >>> d.values()
    [1, 2, 3]
    >>> d['cuatro'] = 4
    >>> d.keys()
    ['uno', 'dos', 'tres', 'cuatro']
    >>> d.items()
    [('uno', 1), ('dos', 2), ('tres', 3), ('cuatro', 4)]
    >>> d.values()
    [1, 2, 3, 4]
    >>> list(d)
    ['uno', 'dos', 'tres', 'cuatro']
    >>> list(d.iterkeys()) == d.keys()
    True
    >>> list(d.itervalues()) == d.values()
    True
    >>> list(d.iteritems()) == d.items()
    True

    """

    def __init__(self, data=()):
        if isinstance(data, GeneratorType):
            data = list(data)
        if isinstance(data, Sequence):
            self.order = [key for key, value in data]
        else:
            self.order = data.keys()
        super(OrderedCaseInsensitiveDict, self).__init__(data)

    def __setitem__(self, key, value):
        if key not in self:
            self.order.append(key)
        super(OrderedCaseInsensitiveDict, self).__setitem__(key, value)

    def __delitem__(self, key):
        raise NotImplementedError

    def __iter__(self):
        return iter(self.order)

    def __deepcopy__(self, memo):
        from copy import deepcopy
        return OrderedCaseInsensitiveDict(
            (key, deepcopy(value, memo)) for key, value in self.iteritems()
        )

    def iterkeys(self):
        return iter(self.order)

    def keys(self):
        return self.order

    def itervalues(self):
        for key in self.order:
            yield self[key]

    def values(self):
        return [self[key] for key in self.order]

    def iteritems(self):
        for key in self.order:
            yield key, self[key]

    def items(self):
        return [(key, self[key]) for key in self.order]


class CaseInsensitiveSet(set):
    """A very basic case-insensitive set.

    >>> s = CaseInsensitiveSet()
    >>> len(s)
    0
    >>> 'a' in s
    False

    >>> list(CaseInsensitiveSet(['aaa', 'Aaa', 'AAA']))
    ['aaa']

    >>> s = CaseInsensitiveSet(['Aaa', 'Bbb'])
    >>> len(s)
    2
    >>> 'aaa' in s
    True
    >>> 'Aaa' in s
    True
    >>> 'AAA' in s
    True
    >>> 'bbb' in s
    True
    >>> 'Bbb' in s
    True
    >>> s.add('ccc')
    >>> len(s)
    3
    >>> 'aaa' in s
    True
    >>> 'ccc' in s
    True
    >>> s.remove('AAA')
    >>> len(s)
    2
    >>> 'aaa' in s
    False

    """

    def __init__(self, *args, **kwargs):
        initial_data = set(*args, **kwargs)
        super(CaseInsensitiveSet, self).__init__(item.lower() for item in initial_data)

    def __contains__(self, item):
        return super(CaseInsensitiveSet, self).__contains__(item.lower())

    def add(self, item):
        super(CaseInsensitiveSet, self).add(item.lower())

    def remove(self, item):
        super(CaseInsensitiveSet, self).remove(item.lower())
