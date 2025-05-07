def get(self, key):
    flags = list(dict(self.all_flags()).keys())
    numbers = list(dict(self.all_doubles()).keys())
    strings = list(dict(self.all_strings()).keys())

    if key in flags:
        return self.get_flag(key)

    if key in numbers:
        value = self.get_numbers(key)
        if len(value) == 1:
            return value[0]
        return value

    if key in strings:
        return self.get_string(key)
        
    raise ValueError(f"unknown configuration parameter: {key}")

def set(self, key, val):
    flags = list(dict(self.all_flags()).keys())
    numbers = list(dict(self.all_doubles()).keys())
    strings = list(dict(self.all_strings()).keys())

    if key in flags:
        return self.set_flag(key, val)

    if key in numbers:
        try:
            len(val)
            return self.set_numbers(key, val)
        except TypeError:
            return self.set_number(key, val)

    if key in strings:
        return self.set_string(key, val)
        
    raise ValueError(f"unknown configuration parameter: {key}")

def set_from_dict(self, d: dict):
    """
    Set config options from a python dictionary.
    """
    for key, val in d.items():
        self.set(key, val)
