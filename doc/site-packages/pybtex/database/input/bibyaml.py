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

import yaml
from pybtex.database.input import BaseParser
from pybtex.database import Entry, Person


class Parser(BaseParser):
    name = 'bibyaml'
    aliases = 'yaml',
    suffixes = '.yaml', '.bibyaml'

    def parse_stream(self, stream):
        t = yaml.safe_load(stream)

        entries = ((key, self.process_entry(entry))
                for (key, entry) in t['entries'].iteritems())

        try:
            self.data.add_to_preamble(t['preamble'])
        except KeyError:
            pass

        self.data.add_entries(entries)
        return self.data

    def process_entry(self, entry):
        e = Entry(entry['type']) 
        for (k, v) in entry.iteritems():
            if k in Person.valid_roles:
                for names in v:
                    e.add_person(Person(**names), k)
            elif k == 'type':
                pass
            else:
                e.fields[k] = unicode(v)
        return e
