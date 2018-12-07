MATHJAX_URL=https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML

README.html: README.rst
	rst2html.py --math-output="MathJax ${MATHJAX_URL}" $^ $@
