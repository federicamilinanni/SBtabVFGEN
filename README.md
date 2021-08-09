<!-- badges: start -->
[![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/a-kramer/SBtabVFGEN/master?urlpath=rstudio)
<!-- badges: end -->

# SBtabVFGEN

Convert a Model written in [SBtab](https://www.sbtab.net/), saved as a
series of `tsv` files or alternatively an [Open Document
Spreadsheet](https://www.documentfoundation.org/) `ods` to a
[VFGEN](https://warrenweckesser.github.io/vfgen/) vector field file
`vf`.

## Install

Using `remotes`:
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

## Shell Interface via `Rscript`

An additional R [script](./R/sbtab_to_vfgen) can be called from the commandline; it uses `sbtab_to_vfgen.R` internally:
```bash
$ ./sbtab_to_vfgen *.tsv
```
This will work if the `.tsv` files have acceptable SBtab content.

## Other Output Formats 

As a by-product, [sbtab_to_vfgen()](./R/sbtab_to_vfgen.R) also produces
a `mod` file intended as a starting point for use in
[neuron](https://neuron.yale.edu/neuron/). This is not the primary
purpose of this function.

If _libSBML_ is installed and R bindings available, an attempt will be
made to produce an SBML file (see further below).

The type of _systems biology_ models that we have in mind are
(plain) ordinary differential equations (ODEs). 

_NEURON_ comprehensively models biochemistry and electrophysiology
(e.g. membrane action potentials). These two different simulations
have to be coupled (neuron does that). This aspect is missing from the
SBtab files we typically use, so the produced `mod` file is really
only an initial point; the user has to change and adapt the file to
the intended purpose and make it work inside of neuron. The user must
also be aware of NEURONs units and make the necessary unit
conversions at some point.

There is no automatic conversion of units (yet).


## VFGEN Output

[VFGEN](https://github.com/WarrenWeckesser/vfgen) is a very useful
tool that reformats an ODE (given in vfgen's `xml` format) and
convert it into various programming languages, including right hand
side functions for two ODE solvers in `C`:
[gsl](https://www.gnu.org/software/gsl/doc/html/ode-initval.html) and
[cvode](https://computing.llnl.gov/projects/sundials/cvode). While
doing so, it uses [GiNaC](https://ginac.de/) to calculate the model's
Jacobian analytically (among other things). 

The R script [sbtab_to_vfgen.R](./sbtab_to_vfgen.R) converts an SBtab model to vfgen's
`xml` format, among others.

Biological models don't necessarily map uniquely onto ODE models, a
compound can be a state variable or an algebraic assignment or a
constant, this has to be inferred a bit from the SBtab files.

## SBtab

[SBtab](https://www.sbtab.net/) is a tabular format for biochemical
models (as in Systems Biology). It is ~perhaps~ easier to understand
than `sbml` and can be parsed/worked via shell scripts (e.g. with line
oriented tools such as `sed` and `awk`) due to its tabular nature
(stored as e.g. tab separated value text files).

In addition to the official upstream documentation, we have summarised the SBtab entries that this script can use in [sbtab.md](./docs/sbtab.md) 

## Open Document Format, Gnumeric and Spreadsheets in General

Even though `.tsv` files are more fundamental types (have least amount
of prerequisites), `.ods` files keep all the sheets in one file and
are slightly more convenient. [Gnumeric](http://www.gnumeric.org/) is
a spreadsheet software that handles both `ods` and tsv files fairly
well (it also has its own format `.gnumeric`)

Conversion between spreadsheet formats like `.ods`, `.gnumeric` and
`.tsv` files is very convenient using `ssconvert`, a part of
[gnumeric](http://www.gnumeric.org/). The shell scripts [ods_to_tsv.sh](./ods_to_tsv.sh) and [tsv_to_ods.sh](./tsv_to_ods.sh) in this
repository are an example of `ssconvert` usage.

An SBtab document can be imported from an open document spreadsheet (`.ods`) directly using the
[readODS](https://cran.r-project.org/web/packages/readODS/index.html)
package:
```R
library("SBtabVFGEN")
model.ods <- "examplemodel.ods" 
if (file.exists(model.ods){
 model.sbtab <- sbtab_from_ods(model.ods)
 sbtab_to_vfgen(model.sbtab)
}
```
The result is written to several files (`.vf`,`.mod`, and `.xml`). Some other results with additional information are also created. 

In either case, whether TSV or ODS was used, `model.sbtab` will be a
list of `data.frame`s.

Other spreadsheet programs such as _google spreadsheets_ and _libre
office_ export to `tsv` _one sheet at a time_ (with no easy
workarounds) and lack an option to export _N_ sheets into _N_
files. Gnumeric's `ssconvert` command does.

### Text files, Unicode and ASCII

Sometimes, spreadsheet software introduces Unicode characters such as
`−` ('MINUS SIGN' U+2212, in html this is «\&minus;») into the
document. They should be replaced by the ascii character `-`
('HYPHEN-MINUS' U+002D) after export to text-based formats. And
similarly for other Unicode characters, unless they appear in comments
(or generally unparsed content). _Minus_ and _Hyphen_ can look quite
similar: `−-` (depending on the chosedn font), but the ascii _hyphen_ is the
character that programming languages understand as _subtraction_.

You can check your files for non-ascii characters like this:
```bash
grep -P '[^[:ascii:]]' *.tsv
#OR
grep -n '[^a-z_A-Z[:digit:][:punct:][:space:]]' *.tsv
# automatic minus to hyphen replacement:
sed -i 's/−/-/g' *.tsv
```

The second option avoids the `-P` switch.

Neither `grep` nor `egrep` define the `[:ascii:]` character class
without the `-P` option for _perl regular expressions_.

Of course you can use [perl](https://www.perl.org/) directly, or
anything else that has regular expressions.

# Systems Biology Markup Language (SBML level 2 version 4)

The program `sbtab_to_vfgen.R` also produces an `.xml` file in the _Systems Biology Markup Language_ (SBML).
This is only done, if _libsbml_ is installed with `R` bindings, like this:

```bash
$ R CMD INSTALL libSBML_5.18.0.tar.gz
```
If this check: `if (require(libSBML))` succeeds, then the scripts attempts to make an sbml file. SBML is a format that has units, and the units defined in SBtab are forwarded to SBML. The formats are very different with regard to unit handling and math generally. The method we use to parse human readble text units is described in [units.md](./units.md).

[Here](./docs/libsbml.md) is a small (incomplete) list of libsbml functions
in R (that we used).

## Some Remarks on LIBSBML

The libSBML R bindings are not documented yet and there is some
guesswork involved (on our side). The code in the repository is not
yet tested very thoroughly, so somewhat cryptic errors may occur. Most
recently we have used the `.tsv` form of sbml, not direct `.ods`
input. Content such as comments, colors, may make the parsing more
difficult. TSV files don't have such _extra_ content, so they have
been more easy to process.

