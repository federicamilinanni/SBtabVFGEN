# SBtabVFGEN

Convert a Model written in [SBtab](https://www.sbtab.net/), saved as an [Open Document Spreadsheet](https://www.documentfoundation.org/) (ods) to a [VFGEN](http://www.warrenweckesser.net/vfgen/) vector field file (vf)

## VFGEN

This is a very useful tool that can reformat an ordinary differential
equation ODE (given in a custom xml format) and convert it into
various programming languages, including right hand side functions for
two `C` solvers of ODEs. While doing so, it uses
[GiNaC](https://ginac.de/) to calculate the model's Jacobian. The R
script `sbtab_to_vfgen` converts an SBtab model to vfgen's xml format.

## SBtab

This is a tabular format for biochemical models (as in
Systems Biology). It is easier to understand than sbml and can be
parsed more easily due to its tabular nature (e.g. csv files).

## Open Document Format

Even though `.[tc]sv` files are perhaps most fundamental, ods files
keep all the sheets in one file and are slightly more convenient. This
script expects the file to be read using the
[readODS](https://cran.r-project.org/web/packages/readODS/index.html)
package.
