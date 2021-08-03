from pybtex.database import parse_file, BibliographyData
import io
import sys
import re

import pybtex.backends.markdown
backend = pybtex.backends.markdown.Backend()

import pybtex.style.formatting.plain
style = pybtex.style.formatting.plain.Style()

def format_list(bib_data):
    formatted = style.format_bibliography(bib_data)

    result = io.StringIO()

    for k, e in enumerate(formatted.entries):
        entry = io.StringIO()
        backend.write_to_stream([e], entry)
        text = entry.getvalue()
        # unwrap lines
        text = text.replace("\n", " ")
        # replace "[N] " with "N. @anchor {key}" to make an ordered list with anchors in
        # Markdown
        text = re.sub(r"^\[([0-9]+)\] ", f"\\1. @anchor {e.key} ", text)

        result.write(text)
        # add a line break at the end of the entry
        result.write("\n")

    return result.getvalue()

header = """References {#references}
==========

@par Notes
This large list collects all references which the PISM authors have found
convenient.  There is no claim that all of these references get direct use,
or even mention, in the PISM project files.

-  -  -  -  -
"""

if __name__ == "__main__":
    input_file = sys.argv[1]

    print(header)
    print(format_list(parse_file(input_file)))
