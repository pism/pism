"""Doctree backend for pybtex."""

import docutils.nodes

from pybtex.backends import BaseBackend
import pybtex.richtext

class Backend(BaseBackend):
    name = 'doctree'

    symbols = {
        'ndash': docutils.nodes.inline(u'\u2013', u'\u2013'),
        'newblock': docutils.nodes.inline(u' ', u' '),
        'nbsp': docutils.nodes.inline(u'\u00a0', u'\u00a0')
    }
    tags = {
         'emph': docutils.nodes.emphasis,
    }

    def format_text(self, text):
        return docutils.nodes.inline(text, text)

    def format_tag(self, tag_name, text):
        tag = self.tags[tag_name]
        if isinstance(text, basestring):
            return tag(text, text)
        else:
            # must be a docutils node
            node = tag('', '')
            node.children.append(text)
            return node

    def format_href(self, url, text):
        if isinstance(url, basestring):
            refuri = url
        else: # isinstance(url, pybtex.richtext.Text)
            refuri = url.plaintext()
        node = docutils.nodes.reference(refuri=refuri)
        node += text
        return node

    def write_entry(self, key, label, text):
        raise NotImplementedError("use Backend.citation() instead")

    def render_sequence(self, text):
        """Return backend-dependent representation of sequence *text*
        of rendered Text objects.
        """
        if len(text) != 1:
            node = docutils.nodes.inline('', '')
            node += text
            return node
        else:
            return text[0]

    def citation(self, entry, document, use_key_as_label=True):
        """Return citation node, with key as name, label as first
        child, and rendered text as second child. The citation is
        expected to be inserted into *document* prior to any docutils
        transforms.
        """
        # see docutils.parsers.rst.states.Body.citation()
        if use_key_as_label:
            label = entry.key
        else:
            label = entry.label
        name = docutils.nodes.fully_normalize_name(entry.key)
        text = entry.text.render(self)
        citation = docutils.nodes.citation()
        citation['names'].append(name)
        citation += docutils.nodes.label('', label)
        citation += text
        document.note_citation(citation)
        document.note_explicit_target(citation, citation)
        return citation

    def citation_reference(self, entry, document, use_key_as_label=True):
        """Return citation_reference node to the given citation. The
        citation_reference is expected to be inserted into *document*
        prior to any docutils transforms.
        """
        # see docutils.parsers.rst.states.Body.footnote_reference()
        if use_key_as_label:
            label = entry.key
        else:
            label = entry.label
        refname = docutils.nodes.fully_normalize_name(entry.key)
        refnode = docutils.nodes.citation_reference(
            '[%s]_' % label, refname=refname)
        refnode += docutils.nodes.Text(label)
        document.note_citation_ref(refnode)
        return refnode
