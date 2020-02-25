from docutils.parsers.rst import Directive, directives
from docutils import nodes
from docutils.core import publish_doctree
from sphinx.errors import SphinxError
import os
import netCDF4

# We store this using a global variable because this much data seems
# to break Sphinx if we store it in an env attribute. (For small
# numbers of parameters config_role an access env via
# inliner.document.settings.env, for larger numbers it cannot.)
pism_data = {}

def make_id(parameter):
    return "config-" + parameter

def format_value(value, T):
    if T == "string" and len(value) == 0:
        return nodes.emphasis("", "no default")

    if T in ["string", "keyword" ]:
        return nodes.literal(value, value.replace(",", ", "))

    if T in ["number", "integer"]:
        return nodes.Text("{:g}".format(value))

    if T == "flag":
        return nodes.Text(str(value))

class config(nodes.literal):
    pass

class softhyphen(nodes.Element):
    pass

def visit_softhyphen_html(self, node):
    self.body.append('&shy;')

def depart_softhyphen_html(self, node):
    pass

def visit_softhyphen_latex(self, node):
    self.body.append(r"\-")

def depart_softhyphen_latex(self, node):
    pass

def config_role(typ, rawtext, text, lineno, inliner, options={}, content=[]):
    """Parse :config:`parameter` roles and create `config` nodes resolved
    in resolve_config_links()."""

    # Warn about misspelled parameters names.
    if text not in pism_data:
        msg = inliner.reporter.error('{}: invalid parameter name'.format(rawtext),
                                     line=lineno)
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]

    return [config("", text)], []

class ParameterList(Directive):
    has_content = False
    required_arguments = 0

    option_spec = {"prefix": directives.unchanged}

    def list_entry(self, name, data):
        "Build an entry for the list of parameters."

        para = nodes.paragraph()
        para += nodes.target('', '', ids=[make_id(name)])
        para += nodes.literal("", name)
        para += nodes.Text(" ({})".format(data["type"]))

        doc, _ = self.state.inline_text(data["doc"], self.lineno)

        fl = nodes.field_list()

        f = nodes.field()
        value = nodes.paragraph()
        value += format_value(data["value"], data["type"])
        if "units" in data:
            value += nodes.emphasis("", " ({})".format(data["units"]))
        f += [nodes.field_name("", "Value"),
              nodes.field_body("", value)]
        fl += f

        if "choices" in data:
            choices = format_value(data["choices"], "keyword")
            f = nodes.field()
            f += [nodes.field_name("", "Choices"),
                  nodes.field_body("", nodes.paragraph("", "", choices))]
            fl += f

        if "option" in data:
            option = format_value("-" + data["option"], "keyword")
            f = nodes.field()
            f += [nodes.field_name("", "Option"),
                  nodes.field_body("", nodes.paragraph("", "", option))]
            fl += f

        f = nodes.field()
        f += [nodes.field_name("", "Description"),
              nodes.field_body("", nodes.paragraph("", "", *doc))]
        fl += f

        para += fl

        return para

    def compact_list_entry(self, name, data):
        "Build an entry for the compact list of parameters."

        doc, _ = self.state.inline_text(data["doc"], self.lineno)

        para = nodes.paragraph()
        para += config("", name)
        para += nodes.Text(" ")
        para += doc
        para += nodes.Text(" (Default: ")
        para += format_value(data["value"], data["type"])
        if "units" in data:
            para += nodes.emphasis("", " ({})".format(data["units"]))

        if "choices" in data:
            para += nodes.Text(", ")
            para += nodes.Text("choices: ")
            choices = format_value(data["choices"], "keyword")
            para += choices

        if "option" in data:
            para += nodes.Text(", ")
            option = format_value("-" + data["option"], "keyword")
            para += nodes.Text("option: ")
            para += option

        para += nodes.Text(".)")

        return para

    def run(self):
        env = self.state.document.settings.env

        # make sure the current document gets re-built if the config file changes
        env.note_dependency(env.pism_parameters["filename"])

        if "prefix" in self.options:
            prefix = self.options["prefix"]
        else:
            prefix = ""

            # Store the docname corresponding to this parameter list.
            # It is used in resolve_config_links() to build URIs.
            env.pism_parameters["docname"] = env.docname

        full_list = prefix == ""

        parameter_list = nodes.enumerated_list()

        parameters_found = False
        for name in sorted(pism_data.keys()):

            # skip parameters that don't have the desired prefix
            if not name.startswith(prefix):
                continue
            else:
                parameters_found = True

            item = nodes.list_item()

            if full_list:
                # only the full parameter list items become targets
                item += self.list_entry(name, pism_data[name])
            else:
                item += self.compact_list_entry(name, pism_data[name])

            parameter_list += item

        if not parameters_found:
            msg = 'Error in a "{}" directive: no parameters with prefix "{}".'.format(self.name, prefix)
            text_error = nodes.error()
            text_error += nodes.Text(msg)

            reporter = self.state_machine.reporter
            system_error = reporter.error(msg, nodes.literal('', ''), line=self.lineno)
            return [text_error, system_error]

        return [parameter_list]

def resolve_config_links(app, doctree, fromdocname):
    """Replace config nodes with references to items in the list of
    parameters.
    """
    docname = app.builder.env.pism_parameters["docname"]

    for node in doctree.traverse(config):
        parameter = node.astext()

        reference = nodes.reference('', '')
        reference['refdocname'] = docname
        reference['refuri'] = "{}#{}".format(app.builder.get_relative_uri(fromdocname, docname),
                                             make_id(parameter))

        words = parameter.split(".")
        reference += nodes.literal("", words[0])
        for w in words[1:]:
            reference += [softhyphen(), nodes.literal("", "." + w)]

        node.replace_self(reference)

def init_pism_parameters(app):
    """Read and pre-process the list of PISM's configuration parameters."""
    global pism_data
    env = app.builder.env

    filename = os.path.normpath(app.config.pism_config_file)

    if not hasattr(env, "pism_parameters"):
        # the path to the configuration file (used to tell Sphinx to
        # re-build if it changes)
        env.pism_parameters = {"filename" : filename}

    f = netCDF4.Dataset(filename)
    variable = f.variables["pism_config"]
    data = {k : getattr(variable, k) for k in variable.ncattrs()}

    suffixes = ["choices", "doc", "option", "type", "units"]

    def special(name):
        "Return true if 'name' is a 'special' parameter, false otherwise."
        for suffix in suffixes:
            if name.endswith("_" + suffix) or name == "long_name":
                return True
        return False

    for key, value in data.items():
        if special(key):
            continue

        pism_data[key] = {"value" : value}

        for s in suffixes:
            try:
                pism_data[key][s] = data[key + "_" + s]
            except:
                pass

def check_consistency(app, env):
    if not "docname" in env.pism_parameters:
        raise SphinxError("make sure that this document contains a pism-parameters directive")

def env_purge_doc(app, env, docname):
    if "docname" in env.pism_parameters and docname == env.pism_parameters["docname"]:
        del env.pism_parameters["docname"]

def setup(app):
    app.add_config_value('pism_config_file', "pism_config.json", 'env')

    app.add_node(config)
    app.add_node(softhyphen,
                 html=(visit_softhyphen_html, depart_softhyphen_html),
                 latex=(visit_softhyphen_latex, depart_softhyphen_latex))
    app.add_role('config', config_role)
    app.add_directive('pism-parameters', ParameterList)

    app.connect("builder-inited", init_pism_parameters)
    app.connect("env-purge-doc", env_purge_doc)
    app.connect("env-check-consistency", check_consistency)
    app.connect('doctree-resolved', resolve_config_links)

    return {
        'version': '0.1',
        'parallel_read_safe': False,
        'parallel_write_safe': True,
    }
