.. include:: ../global.txt

.. default-role:: literal

.. _sec-bug-reports:

Submitting bug reports
======================

Please see the |pism-issues-url| to check if someone already found a similar bug. You can
post an issue there, and label it as a bug, if it is new. Alternatively, send a report by
e-mail to |pism-email|.

Please include the following information in **all** bug-reports and questions about
particular PISM's behavior:

- the PISM version (the output of `pismr -version`),
- the **full** command necessary to reproduce the bug,
- the input files used by the run reproducing the bug,
- a description of what PISM does wrong.

When sending a PISM standard output snippet or a compilation log, please attach it as a
text file instead of pasting it into the message body.

It is great if you can

- remove irrelevant command-line options from the command exhibiting the bug,
- find a quick run that exhibits the bug,
- reduce the size of files it requires, and
- provide a definitive way of checking if the bug in question is fixed.
