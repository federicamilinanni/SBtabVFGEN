# SBtabVFGEN

Convert a Model written in [SBtab](https://www.sbtab.net/), saved as a series of `tsv` files or alternatively an [Open Document Spreadsheet](https://www.documentfoundation.org/) `ods` to a [VFGEN](https://warrenweckesser.github.io/vfgen/) vector field file `vf`.

As a byproduct, this script also produces a `mod` file intended as a starting point for use in [neuron](https://neuron.yale.edu/neuron/).

The type of _systems biology_ models that we have in mind are
autonomous ordinary differential equations (ODEs), while _neuron_ comprehensively models
biochemistry and electrophysiology (action potentials etc.). These two different simulations
have to be coupled (in neuron), so the produced `mod` file is really only
an initial point; the user has to change thye file and make it work inside of
neuron.

## VFGEN

[VFGEN](https://github.com/WarrenWeckesser/vfgen) is a very useful
tool that can reformat an ODE (given in its custom xml format) and
convert it into various programming languages, including right hand
side functions for two `C` solvers of ODEs. While doing so, it uses
[GiNaC](https://ginac.de/) to calculate the model's Jacobian. The R
script `sbtab_to_vfgen.R` converts an SBtab model to vfgen's `xml`
format.

Here is an interactive example:
```R
tsv.files <- dir(pattern=".*[.]tsv");
source("sbtab_to_vfgen.R")
SBtabDoc <- sbtab_from_tsv(tsv.files)
Model <- sbtab_to_vfgen(SBtabDoc)
```

Alternatively, the document can be imported from `ods` using the
[readODS](https://cran.r-project.org/web/packages/readODS/index.html)
package:

```R
ods.file <- "examplemodel.ods" 
source("sbtab_to_vfgen.R")
SBtabDoc <- sbtab_from_ods(ods.file)
Model <- sbtab_to_vfgen(SBtabDoc)
```



## SBtab

[SBtab](https://www.sbtab.net/) is a tabular format for biochemical
models (as in Systems Biology). It is ~perhaps~ easier to understand
than `sbml` and can be parsed via shell scripts (e.g. with line
oriented tools such as `sed` and `awk`) due to its tabular nature
(e.g. tsv files).

We use this format mostly because it can be extended to contain
additional information more easily (such as experimental data, and
conditions under which data was measured).

To make a conversion feasible, we decided on a set of columns and
tables (some specified by `TableName` some also by `TableType`
according to the official specification) which have to be present for
the conversion to work.

In contrast to the official documentation, all tables must be kept in
their own `tsv` file, or different sheets of the same `ods` file. Here
is a list with some information in addition to the obvious `!ID` and
`!Name` columns:

| TableName | Column | Values  | Comment |
| --------: | -----: | :-----: | :------ |
| Compound  | !Scale | log, log10, linear | and some variants of these|
|           | !InitialValue | a number | (per unit) in the above scale |
|           | !Unit | the unit of the above number | as it would be in linear scale |
|           | !SteadyState | `TRUE`/`FALSE` | indicates whether this should converge in a simulation (an output will be generated to monitor the convergence of this item)|
|           | others | unused | but may be informative to the user |
| Parameter | !Scale | log, log10, linear | same as above |
|           | !DefaultValue | a number | in above scale, normnalised to the unit of measurement, possibly subject to fitting/sampling |
|           | !Std | a number | standard deviation / uncertainty of this parameter |
|           | !Min / !Max | numbers | respectively, used if !Std is not present |
| Input     | !DefaultValue | number | same as with normal parameters, but these are set to known values during an _experiment_ (these have to be known values) |
|           | !ConservationLaw | `TRUE`/`FALSE` | this indicates whether a parameter is the result of automatic conservation law analysis (this will be reomved later) |
| Output    | !ErrorName | a string | indicates the column in data sheets that hold the measurement error of an observable (this row) |
|           | !ErrorType | not used | this is for the user |
|           | !ProbDist  | a string | the probability distribution of the noise model (currently unused); this is for humans to read |
|           | !Formula   | a math expression | Outputs are assignments that can be compared to the data (right hand side of assignment) |
| Expression | !Formula | a math expression | right hand side of assignment; Expressions are assignments that are caculated repeatedly each time the ODEs right hand side is called (local variables) |
| Experiments | !Type | `Time␣Series` |  indicates that the data is a `t->output` mapping |
|             |       | `Dose␣Response` | data sheet is an input/output curve, i.e. `input->output` mapping |
|             | `>some_id` | a number | sets the input parameters for this experiment |

All tables require a unique `!ID` column (the ID can be seen as a key
for associative arrays aka _dictionaries_ or _hash tables_). The
`!Name` column must be unique as well and the entries should work as
variable names in the language that you plan to convert the model
to. The script in this repository uses the `make.names()` function on
this column (which will make them unique) and also replaces `.` with
`_`.

Many numbers can be given in a specified scale (like `log`), these
numbers will be converted to linear scale when a model file is written
to file. Let a quantity `y` be measured in unit `M` (y is a number
followed by a unit, y/m is just a number), and !Scale be set to
`log10`, then the number you write in the ![Default]Value column is
`z=log10(y/M)`. The script will do the inverse to generate the model and pass the unit on to `.mod` files.

No unit conversion is attempted here.

A long term goal is that IDs can be used everywhere in the rest of the
document to reference the row in question (but it's not guaranteed
yet) and doesn't yet work in the kinetic law / Reaction Formulas yet
(names are used there). In reference headers `>ID`, the ID must be
used.

The experimental measurement data portion of the files is ignored for
model generation.

If an ID is not unique (but is in a different sheet), new entries from
the new sheet override old entries. IDs have to be unique inside the
same sheet.

## Open Document Format, Gnumeric and Spreadsheets in General

Even though `.tsv` files are more fundamental types (have least amount
of prerequisites), `.ods` files keep all the sheets in one file and
are slightly more convenient.

This `R` script can process such files directly, as mentioned, using
the
[readODS](https://cran.r-project.org/web/packages/readODS/index.html)
package, but you can also convert between spreadsheet formats like
`.ods`, `.gnumeric` and `tsv` files using `ssconvert`, a part of
[gnumeric](http://www.gnumeric.org/). The shell scripts in this
repository are an example of ssconvert usage.

Other spreadsheet programs (like google spreadsheets and libreoffice)
usually export to `tsv` _one sheet at a time_ and lack an option to
export all sheets into as many files.
