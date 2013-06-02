#!/usr/bin/env python

# Copyright (c) 2009, 2010, 2011, 2012  Andrey Golovizin
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

import sys
import optparse

from pybtex.__version__ import version

from pybtex import errors
from pybtex.textutils import capfirst, add_period
from pybtex.plugin import find_plugin, enumerate_plugin_names


def check_plugin(option, option_string, value):
    return find_plugin(option.plugin_group, value)


class PybtexOption(optparse.Option):
    ATTRS = optparse.Option.ATTRS + ['plugin_group']
    TYPES = optparse.Option.TYPES + ('load_plugin',)
    TYPE_CHECKER = dict(optparse.Option.TYPE_CHECKER, load_plugin=check_plugin)


BaseHelpFormatter = optparse.IndentedHelpFormatter
class PybtexHelpFormatter(BaseHelpFormatter):
    def get_plugin_choices(self, plugin_group):
        return ', '.join(sorted(enumerate_plugin_names(plugin_group)))

    def expand_default(self, option):
        result = BaseHelpFormatter.expand_default(self, option)
        if option.plugin_group:
            plugin_choices = self.get_plugin_choices(option.plugin_group)
            result = result.replace('%plugin_choices', plugin_choices)
        return result


make_option = PybtexOption


class CommandLine(object):
    options = ()
    option_defaults = None
    legacy_options = ()
    prog = None
    args = None
    description = ''
    num_args = 0

    def __init__(self):
        self.opt_parser = self.make_option_parser()

    def __call__(self):
        from pybtex.exceptions import PybtexError
        import pybtex.io
        try:
            self.main()
        except PybtexError, error:
            errors.print_error(error)
            sys.exit(1)

    def make_option_parser(self):
        opt_parser = optparse.OptionParser(
            prog=self.prog,
            option_class=PybtexOption,
            formatter=PybtexHelpFormatter(),
            usage='%prog ' + self.args,
            description=capfirst(add_period(self.description)),
            version='%%prog-%s' % version
        )
        for option_group, option_list in self.options:
            if option_group is None:
                container = opt_parser
            else:
                container = opt_parser.add_option_group(option_group)
            for option in option_list:
                container.add_option(option)

        if self.option_defaults:
            opt_parser.set_defaults(**self.option_defaults)

        return opt_parser

    def run(self, options, args):
        raise NotImplementedError

    def recognize_legacy_optons(self, args):
        """Grok some legacy long options starting with a single `-'."""
        return [
            '-' + arg if arg.split('=', 1)[0] in self.legacy_options else arg
            for arg in args
        ]

    def main(self):
        args = self.recognize_legacy_optons(sys.argv[1:])
        options, args = self.opt_parser.parse_args(args)
        if len(args) != self.num_args:
            self.opt_parser.print_help()
            sys.exit(1)

        self.run(options, args)
        sys.exit(errors.error_code)
