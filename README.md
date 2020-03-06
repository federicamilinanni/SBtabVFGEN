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
their own `tsv` file, or different sheets of the same `ods` file. The
following Sections have more information on specific tables and their
required columns (in addition to the obvious `!ID` and `!Name`
columns).

All tables require a unique `!ID` column (the ID can be seen as a key
for associative arrays aka _dictionaries_ or _hash tables_). The
`!Name` column must be unique as well and the entries should work as
variable names in the language that you plan to convert the model
to. The script in this repository uses the `make.names()` function on
this column (which will make them unique, but break reactions) and
also replaces `.` with `_` in all names (dots are often illegal in
variable names, not in `R` though).

Many numbers can be given in a specified scale (like `log`), these
numbers will be converted to linear scale when a model file is written
to file. Let a quantity `y` be measured in unit `M` (y is a number
followed by a unit, y/m is just a number), and !Scale be set to
`log10`, then the number you write in the ![Default]Value column is
`z=log10(y/M)`. The script will do the inverse to generate the model and pass the unit on to `.mod` files.

No unit conversion is attempted here.

The names of all tables must be unique.

### Compound

This table defines the compounds that are supposed to be modeled bu
state variables and are subject to change by the reactions in the
systems. Currently, compounds that are unchangeable (fixed by external
conditions) should not be here.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !Scale | log, log10, linear | and some variants of these|
| !InitialValue | a number | (per unit) in the above scale |
| !Unit | the unit of the above number | as it would be in linear scale |
| !SteadyState | `TRUE` | this compound should reach a steady state in at least onw scenario and you want to know whether this happened |
|              |`FALSE` | it is not important whether or not this compound reaches steady state|

The conversion script will make a file called
`[…]SuggestedOutput.tsv`, it will have lines that can be used to check
whether a compound has reached steady state (or not), this is done for
each compound that has `!SteadyState` marked as `TRUE`. If that output
is close to `0`, then steady state was reached (it's the sum of all
fluxes for the compound in question).

Others columns are unused but may be informative to the user, or others.

### Parameters

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !Scale | `log`, `log10`, `linear` | some aliases of these are possible (such as `base-10 logarithm`)|
| !DefaultValue | a number | in above scale, normnalised to the unit of measurement, possibly subject to fitting/sampling |
| !Std | a number | standard deviation / uncertainty of this parameter |
| !Min / !Max | numbers | respectively, used if !Std is not present |

The columns `!Std` and `!Min/Max` are only used in
sampling/optimisation, the model conversion is unaffected by them, the
DefaultValue is passed on to the model files (if there is a place to
put them).

### Reactions

The column `!ReactionFormula` determines the stoichiometry of the
model, the `!KineticLaw` column determines the flux of the given
reaction. Both er required and are standard columns in SBtab.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !KineticLaw | e.g. `kf*A*B-kr*AB` | the flux, as a math expression |
| !ReactionFormula | e.g. `A+B<=>AB` | so, `AB` will increase and both `A` and `B` will decrease by this reaction whenever the flux is positive |

Since the kinetic law dtermines the reversibility of the reaction, the
column `!IsReversible` is not necessary, but if you determine the
kinetics based on the law of mass action it may be important for you
to have that column as a reminder (for when you are auto generating
the `!KineticLaw` column, which this script doesn't do).

### Input

The input parameters to the model that distinguish different
experiments. These quantities are known and can be influenced by the
people who are performaing an experiment (or rather the real
counterparts of these quantities can be influenced). These play the
roles of (additional) parameters, but a different kind of parameter
than in the Parameter table. Experiments are supposed to have the same
parameters of the normal kind and different parameters of the input
kind.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !DefaultValue | number | same as with normal parameters, but these are set to known values during an _experiment_ (these have to be known values) |
| !ConservationLaw | `TRUE`| the parameter is the result of an earlier run of the conversion script|
||`FALSE` | this parameter is unrelated to conservation laws |

The `!ConservationLaw` column will eventually be obsolete and the
numbers in that column will be determined by the _experiment-specific_
initial conditions. Currently there is only one _initial condition_ vector for
_all_ experiments.

At this moment, a conservation law parameter will be ignored whenever
the script runs (again), because the script creates those itself (every time it runs).

### Output

The outputs are _observable quantities_ of this system; what is and
isn't an output depends on what you can measure (or have knowlege
about). Outputs are ususally converted to functions in the target
language. _Experimental Data_ and _Outputs_ are intimately related as
the outputs are the model's equivalent of the data and in some way
those can be compared to one another.

It is possible to include the measured data in other sheets, that data should be
stored together with an estimate of the measurement noise
levels. Regardless of the nature of the noise and underlying
distributions we use `!ErrorName` to indicate which column (the one
that has this name) is storing information about this measurement
error.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !ErrorName | a string | indicates the column in data sheets that hold the measurement error of an observable |
| !ErrorType | not used | this is for the user |
| !ProbDist  | a string | the probability distribution of the noise model (currently unused); this is for humans to read |
| !Formula   | a math expression | the right hand side of the assignment |

Many data columns may share the same Error column. This is useful if
you have only a very rough estimate of the noise anyway, and the
outputs are in the same number range so using the same standard
deviation (etc.) for all data points seems good enough.

The data sheets are not used by this script, but they are used [here](https://github.com/a-kramer/mcmc_clib).

### Expression

These will be local variables of the model. Expressions are
assignments that are caculated repeatedly each time the ODEs right
hand side is called (before fluxes are calculated, fluxes are otherwise also a
type of _expression_).

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !Formula | a math expression | right hand side of assignment|

Currently, there is no way to formally write something like _initial
assignments_, as those don't go into the model files, these
assignments have to happen before you call the model solver, so those
are up to you use of the model in the language you had in mind.

The only type of initial assignment that can be specified is the state
variables' inital conditions.

### Experiments

This table holds the mapping between input parameters and data
sheets. It determines the conditions under which a data set should be
replicated using the model.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !Type | `Time␣Series` |  indicates that the data is a `t->output` mapping |
|       | `Dose␣Response` | data sheet is an input/output curve, i.e. `input->output` mapping |
| `>some_id` | a number | sets the input parameters for this experiment |
|!Event| a table name | the name of an event table that holds time instantaneous model state changes |

This is not used by this converter, but useful for parameter fitting
and interpretation of _the input_ and _output_ concepts.

### Remarks about `!ID`s and `!Name`s

A long term goal is that IDs can be used everywhere in the rest of the
document to reference the row in question (but it's not guaranteed
yet) and doesn't yet work in the kinetic law / Reaction Formulas
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
