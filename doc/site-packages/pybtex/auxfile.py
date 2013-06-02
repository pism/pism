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

"""parse LaTeX aux file
"""

from __future__ import with_statement

import re

from pybtex.exceptions import PybtexError
import pybtex.io


class AuxDataError(PybtexError):
    pass


class AuxData:
    command_re = re.compile(r'\\(citation|bibdata|bibstyle|@input){(.*)}')
    def __init__(self, encoding):
        self.filename = None
        self.encoding = encoding
        self.style = None
        self.data = None
        self.citations = []

    def handle_citation(self, s):
        for c in s.split(','):
            if not c in self.citations:
                self.citations.append(c)

    def handle_bibstyle(self, style):
        if self.style is not None:
            raise AuxDataError(r'illegal, another \bibstyle command in %s' % self.filename)
        self.style = style

    def handle_bibdata(self, bibdata):
        if self.data is not None:
            raise AuxDataError(r'illegal, another \bibdata command in %s' % self.filename)
        self.data = bibdata.split(',')

    def handle_input(self, filename):
        self.parse_file(filename)

    def handle(self, command, value):
        action = getattr(self, 'handle_%s' % command.lstrip('@'))
        action(value)

    def parse_file(self, filename):
        previous_filename = self.filename
        self.filename = filename

        with pybtex.io.open_unicode(filename, encoding=self.encoding) as f:
            s = f.read()
        for command, value in self.command_re.findall(s):
            self.handle(command, value)

        self.filename = previous_filename


def parse_file(filename, encoding):
    """Parse a file and return an AuxData object."""

    data = AuxData(encoding)
    data.parse_file(filename)
    if data.data is None:
        raise AuxDataError(r'found no \bibdata command in %s' % filename)
    if data.style is None:
        raise AuxDataError(r'found no \bibstyle command in %s' % filename)
    return data
