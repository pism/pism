import os
import pkgutil
import posixpath
from contextlib import contextmanager
from shutil import rmtree
from tempfile import mkdtemp


from pybtex import io
from pybtex import errors
from pybtex import bibtex
from pybtex.tests import diff


@contextmanager
def cd_tempdir():
    current_workdir = os.getcwd()
    tempdir = mkdtemp(prefix='pybtex_test_')
    os.chdir(tempdir)
    try:
        yield tempdir
    finally:
        os.chdir(current_workdir)
        rmtree(tempdir)


def copy_resource(package, resource):
    filename = posixpath.basename(resource)
    data = pkgutil.get_data(package, resource).decode(io.get_default_encoding())
    with io.open_unicode(filename, 'w') as data_file:
        data_file.write(data)


def copy_files(bib_name, bst_name):
    copy_resource('pybtex.tests.data', bib_name + '.bib')
    copy_resource('pybtex.tests.data', bst_name + '.bst')


def write_aux(aux_name, bib_name, bst_name):
    with io.open_unicode(aux_name, 'w') as aux_file:
        aux_file.write(u'\\citation{*}\n')
        aux_file.write(u'\\bibstyle{{{0}}}\n'.format(bst_name))
        aux_file.write(u'\\bibdata{{{0}}}\n'.format(bib_name))


def check_make_bibliography(bib_name, bst_name):
    with cd_tempdir() as tempdir:
        copy_files(bib_name, bst_name)
        write_aux('test.aux', bib_name, bst_name)
        with errors.capture() as stderr:  # FIXME check error messages
            bibtex.make_bibliography('test.aux')
        with io.open_unicode('test.bbl', 'r') as result_file:
            result = result_file.read()
        correct_result_name = '{0}_{1}.bbl'.format(bib_name, bst_name)
        correct_result = pkgutil.get_data('pybtex.tests.data', correct_result_name).decode(io.get_default_encoding())
        assert result == correct_result, diff(correct_result, result)


def test_bibtex_engine():
    for bib_name, bst_name in [
        ('xampl', 'unsrt'),
        ('xampl', 'plain'),
        ('cyrillic', 'unsrt'),
    ]:
        yield check_make_bibliography, bib_name, bst_name
