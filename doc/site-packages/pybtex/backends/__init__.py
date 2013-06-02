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

import pybtex.io
from pybtex.plugin import Plugin


available_plugins = ('latex', 'html', 'plaintext', 'doctree')


class BaseBackend(Plugin):
    default_plugin = 'latex'

    def __init__(self, encoding=None):
        self.encoding = encoding

    def write_prologue(self):
        pass

    def write_epilogue(self):
        pass

    def format_text(self, text):
        return text

    def format_tag(self, tag_name, text):
        """Format a tag with some text inside.

        Text is already formatted with format_text."""

        raise NotImplementedError

    def format_href(self, url, text):
        """Format a hyperlink with some text inside.

        Text is already formatted with format_text."""

        raise NotImplementedError

    def render_sequence(self, text):
        """Render a sequence of rendered text objects."""
        return "".join(text)

    def write_entry(self, label, key, text):
        raise NotImplementedError

    def write_to_file(self, formatted_entries, filename):
        with pybtex.io.open_unicode(filename, "w", self.encoding) as stream:
            self.write_to_stream(formatted_entries, stream)

    def write_to_stream(self, formatted_bibliography, stream):
        self.output = stream.write
        self.formatted_bibliography = formatted_bibliography

        self.write_prologue()
        for entry in formatted_bibliography:
            self.write_entry(entry.key, entry.label, entry.text.render(self))
        self.write_epilogue()

    def write_bibliography(self, formatted_bibliography, filename):
        import warnings
        warnings.warn(
            'write_bibliography() is deprecated, use write_to_file() insted.',
            DeprecationWarning,
        )
        self.write_to_file(formatted_bibliography, filename)
