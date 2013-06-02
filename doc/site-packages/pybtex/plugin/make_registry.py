#!/usr/bin/env python

import os
import imp
import pkgutil
import json


plugin_types = (
    ('pybtex.database.input', 'Parser'),
    ('pybtex.database.output', 'Writer'),
    ('pybtex.backends', 'Backend'),
    ('pybtex.style.labels', 'LabelStyle'),
    ('pybtex.style.names', 'NameStyle'),
    ('pybtex.style.sorting', 'SortingStyle'),
    ('pybtex.style.formatting', 'Style'),
)


def import_from(module_name, obj_name):
    module = __import__(module_name, globals(), locals(), [obj_name])
    return module, getattr(module, obj_name)


def iter_plugins(base_module, class_name):
    plugin_dir = os.path.dirname(base_module.__file__)
    for loader, name, is_pkg in pkgutil.iter_modules([plugin_dir]):
        plugin_module = imp.load_module(name, *imp.find_module(name, [plugin_dir]))
        yield getattr(plugin_module, class_name)


def get_plugin_group_info(plugin_group, class_name):
    base_module, base_plugin = import_from(plugin_group, 'Base' + class_name)
    info = {
        'class_name': class_name,
        'default_plugin': base_plugin.default_plugin,
        'plugins': [],
        'aliases': {},
        'suffixes': {},
    }
    for plugin in iter_plugins(base_module, class_name):
        info['plugins'].append(plugin.name)
        info['aliases'].update(dict.fromkeys(plugin.aliases, plugin.name))
        info['suffixes'].update(dict.fromkeys(plugin.suffixes, plugin.name))
    return info


def get_plugin_registry():
    return dict(
        (plugin_group, get_plugin_group_info(plugin_group, class_name))
        for plugin_group, class_name in plugin_types
    )


def write_plugin_registry(filename):
    with open(filename, 'wb') as registry_file:
        registry_file.write('plugin_registry = ')
        json.dump(get_plugin_registry(), registry_file, indent=4)
        registry_file.write('\n')


if __name__ == '__main__':
    write_plugin_registry('registry.py')
