# vim:fileencoding=utf-8

from pybtex.database import Entry, Person
from pybtex.database import BibliographyData

reference_data = BibliographyData(
    entries=[
        ('ruckenstein-diffusion', Entry('article',
            fields={
                'language': u'english',
                'title': u'Predicting the Diffusion Coefficient in Supercritical Fluids',
                'journal': u'Ind. Eng. Chem. Res.',
                'volume': u'36',
                'year': u'1997',
                'pages': u'888-895',
            },
            persons={'author': [Person(u'Liu, Hongquin'), Person(u'Ruckenstein, Eli')]},
        )),
        ('test-booklet', Entry('booklet',
            fields={
                'language': u'english',
                'title': u'Just a booklet',
                'year': u'2006',
                'month': u'January',
                'address': u'Moscow',
                'howpublished': u'Published by Foo',
            },
            persons={'author': [Person(u'de Last, Jr., First Middle')]}
        )),
        ('test-inbook', Entry('inbook',
            fields={
                'publisher': u'Some Publisher',
                'language': u'english',
                'title': u'Some Title',
                'series': u'Some series',
                'booktitle': u'Some Good Book',
                'number': u'3',
                'edition': u'Second',
                'year': u'1933',
                'pages': u'44--59',
            },
            persons={'author': [Person(u'Jackson, Peter')]}
        )),
        ('viktorov-metodoj', Entry('book',
            fields={
                'publisher': u'Л.: <<Химия>>',
                'year': u'1977',
                'language': u'russian',
                'title': u'Методы вычисления физико-химических величин и прикладные расчёты',
            },
            persons={'author': [Person(u'Викторов, Михаил Маркович')]}
        )),
    ],
    preamble=['%%% pybtex example file']
)
