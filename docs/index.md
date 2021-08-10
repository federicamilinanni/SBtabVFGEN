# SBtabVFGEN

This is an R package which converts a Model written in
[SBtab](https://www.sbtab.net/), saved as a series of `tsv` files or
alternatively an [Open Document
Spreadsheet](https://www.documentfoundation.org/) `ods` to a
[VFGEN](https://warrenweckesser.github.io/vfgen/) vector field file
`vf`.

## Install

Using `remotes` (R):
```R
remotes::install_github("a-kramer/SBtabVFGEN")
```
You may have to check `.libPaths()` to verify that it includes a path
that you have permission to write to (this is just generally the case,
not just for this package).

## Usage Example

Within an interactive R session called from the folder that contains
the tsv files:
```R
library(SBtabVFGEN)
model.tsv <- dir(pattern=".*[.]tsv$");
model.sbtab <- sbtab_from_tsv(model.tsv)
sbtab_to_vfgen(model.sbtab)
```

# SBtab (Systems Biology Tables)

This format is human readable and writable. We use the standard as we
understand it, but use our own code to process the files, using a
[subset of named columns](./sbtab.md). These files can hold a great
variety of content.

# Systems Biology Markup Language (SBML level 2 version 4)

The program in `sbtab_to_vfgen.R` also produces an `.xml` file in the _Systems Biology Markup Language_ (SBML).
This is only done, if _libsbml_ is installed with `R` bindings, like this:

```bash
$ R CMD INSTALL libSBML_5.18.0.tar.gz
```

If this check: `if (require(libSBML))` succeeds, then the scripts
attempts to make an sbml file. SBML is a format that has _units_.  The
units found in the SBtab document (strings) are forwarded to SBML
(nested, structured xml-elements). The formats are very different with
regard to unit handling and math generally. The method we use to parse the
common (human readble) text units is described in [units.md](./units.md).

There is an official
[guide](http://sbml.org/Software/libSBML/libSBML_R_Example_Programs)
on [libsbml.org](libsbml.org) that we used to write the SBML output
part of this program. The guide contains some examples.
But, the libSBML R bindings are not documented yet and there is some
guesswork involved (on our side). 

[Here](./docs/libsbml.md) is a small (incomplete) list of libsbml functions
in R that we use to create SBML files.

