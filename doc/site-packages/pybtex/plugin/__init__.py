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

import os
import sys
from itertools import chain

from pybtex.exceptions import PybtexError


class PluginGroupNotFound(PybtexError):
    def __init__(self, group_name):
        message = u'plugin group {group_name} not found'.format(
            group_name=group_name,
        )
        super(PluginGroupNotFound, self).__init__(message)


class PluginNotFound(PybtexError):
    def __init__(self, plugin_group, name=None, filename=None):
        assert plugin_group
        if name:
            message = u'plugin {plugin_group}.{name} not found'.format(
                plugin_group=plugin_group,
                name=name,
            )
        elif filename:
            message = u'cannot determine file type for {filename}'.format(
                plugin_group=plugin_group,
                filename=filename,
            )
        else:
            message = u'' # FIXME

        super(PluginNotFound, self).__init__(message)


class Plugin(object):
    name = None
    aliases = ()
    suffixes = ()
    default_plugin = None

    @classmethod
    def get_default_suffix(cls):
        return cls.suffixes[0]


class PluginLoader(object):
    def find_plugin(plugin_group, name):
        raise NotImplementedError

    def enumerate_plugin_names(self, plugin_group):
        raise NotImplementedError


class BuiltInPluginLoader(PluginLoader):
    def get_group_info(self, plugin_group):
        from pybtex.plugin.registry import plugin_registry
        try:
            return plugin_registry[plugin_group]
        except KeyError:
            raise PluginGroupNotFound(plugin_group)

    def find_plugin(self, plugin_group, name=None, filename=None):
        plugin_group_info = self.get_group_info(plugin_group)
        module_name = None
        if name:
            if name in plugin_group_info['plugins']:
                module_name = name
            elif name in plugin_group_info['aliases']:
                module_name = plugin_group_info['aliases'][name]
            else:
                raise PluginNotFound(plugin_group, name=name)
        elif filename:
            suffix = os.path.splitext(filename)[1]
            if suffix in plugin_group_info['suffixes']:
                module_name = plugin_group_info['suffixes'][suffix]
            else:
                raise PluginNotFound(plugin_group, filename=filename)
        else:
            module_name = plugin_group_info['default_plugin']

        class_name = plugin_group_info['class_name']
        return self.import_plugin(plugin_group, module_name, class_name)

    def import_plugin(self, plugin_group, module_name, class_name):
        module = __import__(str(plugin_group), globals(), locals(), [str(module_name)])
        return getattr(getattr(module, module_name), class_name)

    def enumerate_plugin_names(self, plugin_group):
        plugin_group_info = self.get_group_info(plugin_group)
        return plugin_group_info['plugins']


class EntryPointPluginLoader(PluginLoader):
    def find_plugin(self, plugin_group, name=None, filename=None):
        try:
            import pkg_resources
        except ImportError:
            raise PluginNotFound(plugin_group)

        def load_entry_point(group, name):
            entry_points = pkg_resources.iter_entry_points(group, name)
            try:
                entry_point = entry_points.next()
            except StopIteration:
                raise PluginNotFound(plugin_group, name, filename)
            else:
                return entry_point.load()

        if name:
            return load_entry_point(plugin_group, name)
        elif filename:
            suffix = os.path.splitext(filename)[1]
            return load_entry_point(plugin_group + '.suffixes', suffix)
        else:
            raise PluginNotFound(plugin_group)

    def enumerate_plugin_names(self, plugin_group):
        try:
            import pkg_resources
        except ImportError:
            return
        entry_points = pkg_resources.iter_entry_points(plugin_group)
        return [entry_point.name for entry_point in entry_points]


plugin_loaders = [EntryPointPluginLoader(), BuiltInPluginLoader()]


def find_plugin(plugin_group, name=None, filename=None):
    if isinstance(name, type) and issubclass(name, Plugin):
        plugin = name
        #assert plugin.group_name == plugin_group
        return plugin
    else:
        for loader in plugin_loaders:
            try:
                return loader.find_plugin(plugin_group, name, filename)
            except PluginNotFound:
                continue
        raise PluginNotFound(plugin_group, name, filename)


def enumerate_plugin_names(plugin_group):
    results = (
        loader.enumerate_plugin_names(plugin_group)
        for loader in plugin_loaders
    )
    return chain(*results)
