plugin_registry = {
    "pybtex.database.output": {
        "class_name": "Writer", 
        "suffixes": {
            ".xml": "bibtexml", 
            ".bibtexml": "bibtexml", 
            ".bibyaml": "bibyaml", 
            ".bib": "bibtex", 
            ".yaml": "bibyaml"
        }, 
        "aliases": {
            "yaml": "bibyaml"
        }, 
        "default_plugin": "bibtex", 
        "plugins": [
            "bibtex", 
            "bibtexml", 
            "bibyaml"
        ]
    }, 
    "pybtex.style.formatting": {
        "class_name": "Style", 
        "suffixes": {}, 
        "aliases": {}, 
        "default_plugin": "unsrt", 
        "plugins": [
            "plain", 
            "unsrt",
            "alpha",
            "unsrtalpha",
        ]
    }, 
    "pybtex.style.labels": {
        "class_name": "LabelStyle", 
        "suffixes": {}, 
        "aliases": {}, 
        "default_plugin": "number", 
        "plugins": [
            "number",
            "alpha",
        ]
    }, 
    "pybtex.backends": {
        "class_name": "Backend", 
        "suffixes": {
            ".html": "html", 
            ".txt": "plaintext", 
            ".bbl": "latex"
        }, 
        "aliases": {
            "text": "plaintext"
        }, 
        "default_plugin": "latex", 
        "plugins": [
            "html", 
            "latex", 
            "plaintext"
        ]
    }, 
    "pybtex.database.input": {
        "class_name": "Parser", 
        "suffixes": {
            ".xml": "bibtexml", 
            ".bibtexml": "bibtexml", 
            ".bibyaml": "bibyaml", 
            ".bib": "bibtex", 
            ".yaml": "bibyaml"
        }, 
        "aliases": {
            "yaml": "bibyaml"
        }, 
        "default_plugin": "bibtex", 
        "plugins": [
            "bibtex", 
            "bibtexml", 
            "bibyaml"
        ]
    }, 
    "pybtex.style.names": {
        "class_name": "NameStyle", 
        "suffixes": {}, 
        "aliases": {
            "last_first": "lastfirst"
        }, 
        "default_plugin": "plain", 
        "plugins": [
            "lastfirst", 
            "plain"
        ]
    }, 
    "pybtex.style.sorting": {
        "class_name": "SortingStyle", 
        "suffixes": {}, 
        "aliases": {}, 
        "default_plugin": "none", 
        "plugins": [
            "author_year_title", 
            "none"
        ]
    }
}
