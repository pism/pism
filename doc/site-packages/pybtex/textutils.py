# Copyright (c) 2006, 2007, 2008, 2009, 2010, 2011, 2012  Andrey Golovizin
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

import re

terminators = '.?!'
dash_re = re.compile(r'-')
whitespace_re = re.compile(r'\s+')

def capfirst(s):
    return s[0].upper() + s[1:] if s else s

def is_terminated(s):
    """Return true if s ends with a terminating character.
    """
    return (bool(s) and s[-1] in terminators)

def add_period(s):
    """Add a period to the end of s, if there is none yet.
    """
    if s and not is_terminated(s):
        return s + '.'
    return s

def abbreviate(s):
    """Abbreviate some text.
    Examples:
    abbreviate('Some words') -> "S. w."
    abbreviate('First-Second') -> "F.-S."
    """
    def parts(s):
        start = 0
        length = 0
        for letter in s:
            length += 1
            if not letter.isalpha():
                yield s[start:length], letter
                start += length
                length = 0
        yield s[start:length], ""
    def abbr(part):
        if part[0]:
            if is_terminated(part[1]):
                return part[0][0].upper() + part[1]
            else:
                return part[0][0].upper() + '.'
        else:
            return part[1]
    return ''.join(abbr(part) for part in parts(s))

def normalize_whitespace(string):
    r"""
    Replace every sequence of whitespace characters with a single space.

    >>> print normalize_whitespace('abc')
    abc
    >>> print normalize_whitespace('Abc def.')
    Abc def.
    >>> print normalize_whitespace(' Abc def.')
    Abc def.
    >>> print normalize_whitespace('Abc\ndef.')
    Abc def.
    >>> print normalize_whitespace('Abc\r\ndef.')
    Abc def.
    >>> print normalize_whitespace('Abc    \r\n\tdef.')
    Abc def.
    >>> print normalize_whitespace('   \nAbc\r\ndef.')
    Abc def.
    """

    return whitespace_re.sub(' ', string.strip())
