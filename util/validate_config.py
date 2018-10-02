#!/usr/bin/env python

import json
import jsonschema
import copy

# Shortcuts to save some typing:
integer_t = {"type": "integer"}
number_t = {"type": "number"}
string_t = {"type": "string"}
boolean_t = {"type": "boolean"}

# A generic entry:
entry = {
    "type": "object",
    "properties": {
        "doc": string_t,
        "reference": string_t,
        "comment": string_t
    },
    "required": ["doc", "value", "type"],
    "additionalProperties": False
}


def gen_type(keyword):
    return {"type": "string", "pattern": "^%s$" % keyword}


# numerical configuration parameters
number = copy.deepcopy(entry)
number['properties']['value'] = number_t
number['properties']['units'] = string_t
number['properties']['type'] = gen_type("number")
number['required'].append('units')

# integer configuration parameters
integer = copy.deepcopy(number)
integer['properties']['value'] = integer_t
integer['properties']['type'] = gen_type("integer")

# strings
string = copy.deepcopy(entry)
string['properties']['value'] = string_t
string['properties']['type'] = gen_type("string")

# keywords ("choose one of...")
keyword = copy.deepcopy(string)
keyword['properties']['choices'] = {
    "type": "array",
    "minItems": 2,
    "items": {
        "type": "string"
    }
}

# configuration flags
boolean = copy.deepcopy(entry)
boolean['properties']['type'] = gen_type("boolean")
boolean['properties']['value'] = boolean_t

# assembled schema
config_schema = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "type": "object",
    "definitions": {
        "number": number,
        "integer": integer,
        "string": string,
        "keyword": keyword,
        "boolean": boolean,
    },
    "properties": {
        "doc": string_t,
        "comment": string_t,
        "reference": string_t,
    },
    "required": ["doc"],
    "additionalProperties": {
        "anyOf": [
            {"$ref": "#"},
            {"$ref": "#/definitions/number"},
            {"$ref": "#/definitions/integer"},
            {"$ref": "#/definitions/string"},
            {"$ref": "#/definitions/keyword"},
            {"$ref": "#/definitions/boolean"},
        ]
    }
}


def reload_data():
    f = open("pism_config.json")
    j = json.load(f)
    f.close()

    return j


j = reload_data()


def validate_recursively(data, path=[]):
    try:
        jsonschema.validate(data, config_schema)
    except jsonschema.ValidationError as e:
        if len(e.path) == 0:
            print("Failed while validating %s" % json.dumps(e.instance))
            print("Message:", e.message)
            print("Path:", ".".join(path))
        else:
            validate_recursively(e.instance, path + list(e.path))


validate_recursively(j)
