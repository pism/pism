#!/usr/bin/env python3

import csv
import time

year = time.gmtime(time.time())[0]

funding = {}
with open("funding.csv", "r") as f:
    reader = csv.reader(f, skipinitialspace=True, quoting=csv.QUOTE_ALL)

    funding = {}
    for row in reader:
        start_year, end_year, agency, number, _ = row

        try:
            start_year = int(start_year)
            end_year = int(end_year)
        except:
            continue

        # skip grants for which we don't have a number (yet)
        if number.strip() == "":
            continue

        if start_year <= year and year <= end_year:
            try:
                funding[agency].append(number)
            except:
                funding[agency] = [number]

def join(strings):
    assert len(strings) > 0
    if len(strings) == 1:
        return strings[0]
    elif len(strings) == 2:
        return "{} and {}".format(strings[0], strings[1])
    else:
        return join(["{}, {}".format(strings[0], strings[1]),
                     join(strings[2:])])

grants = []
for k, v in funding.items():
    suffix = "s" if len(v) > 1 else ""

    grants.append("{agency} grant{suffix} {number}".format(agency=k,
                                                           suffix=suffix,
                                                           number=join(v)))

print("""
..
   DO NOT EDIT: This file was automatically generated by running doc/funding.py

   Edit doc/funding.py and doc/funding.csv instead.
""")

print("""
Development of PISM is supported by {grants}.""".format(grants=join(grants)))