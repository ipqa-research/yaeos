project: Yaeos
summary: Yet another Equation of State (library)
project_github: https://github.com/fedebenelli/yaeos
author: Federico Benelli
author_description: PhD student with focus on reservoir PVT simulation.
author_email: federico.benelli@mi.unc.edu.ar
github: https://github.com/fedebenelli
src_dir: ../src
exclude_dir: ../test ../doc
output_dir: ../doc/ford_site
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
page_dir:pages

{!../README.md!}
