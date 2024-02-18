---
project: yaeos
summary: Calculation of thermodynamic properties and phase-equilibria with
         Equation of State. Keeping high performance with the power of Fortran
         and the possiblity of using either automatic differentiation and
         anallytical derivatives.
project_github: https://github.com/fedebenelli/yaeos
author: Federico Benelli
author_description: PhD student with focus on reservoir PVT simulation.
email: federico.benelli@mi.unc.edu.ar
github: https://github.com/fedebenelli
src_dir: src
exclude_dir: test doc
output_dir: doc/ford_site
preprocessor: gfortran -E
display: public
         protected
         private
source: false
proc_internals: true
sort: permission-alpha
docmark_alt: !>
docmark: !
predocmark_alt: *
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
graph: false
license: MPL
page_dir: doc/page
media_dir: doc/media
---

[TOC]

{!README.md!}
