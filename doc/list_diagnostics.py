import json
from argparse import ArgumentParser

parser = ArgumentParser()
parser.description = '''Generate an Org table documenting PISM's diagnostics.'''
parser.add_argument("FILE", nargs=1)
parser.add_argument("--type", dest="output_type", choices=["org"],
                    default="org")
options = parser.parse_args()

with open(options.FILE[0]) as f:
    diagnostics = json.load(f)

def print_diagnostics_org():
    def print_table(title, diagnostics):
        print("""
** {title}

#+caption: {title}
| Diagnostic | NetCDF Variable | Units | Description | CF Standard Name | Comment |
|------------|-----------------|-------|-------------|------------------|---------|""".format(title=title))

        names = diagnostics.keys()
        names.sort()
        for name in names:
            for data in diagnostics[name]:
                var_name, units, long_name, standard_name, comment = data
                if len(units) == 0:
                    units = "not set"
                if len(standard_name) == 0:
                    standard_name = "none"
                print('| ={}= | ={}= | ={}= | {} | ={}= | {} |'.format(name, var_name, units, long_name, standard_name, comment))

    print("#+options: ^:nil")

    print("* List of PISM's diagnostics")

    print_table("Spatial diagnostics", diagnostics["spatial"])
    print_table("Scalar diagnostics",  diagnostics["scalar"])

if options.output_type == "org":
    print_diagnostics_org()
