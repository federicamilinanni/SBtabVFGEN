# SBtab

We use this format mostly because it can be extended to contain
additional information more easily (such as experimental data, and
conditions under which data was measured).

To make a conversion feasible, we decided on a set of columns and
tables (some specified by `TableName` some also by `TableType`
according to the official specification) which have to be present for
the conversion to work.

In contrast to the [official
documentation](https://www.sbtab.net/sbtab/default/documentation.html),
we need all tables to be kept in their own `.tsv` file (not all tables
in one huge tsv file), or different sheets of the same `ods` file (the
official SBtab project uses `.xlsx`). But, we don't use any of the
code from theoriginal [SBtab authors](https://www.sbtab.net/sbtab/default/team.html).

The following Sections have more information on specific tables and
their required columns (in addition to the obvious `!ID` and `!Name`
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
followed by a unit, y/M is just a number), and `!Scale` be set to
`log10`, then the number you write in the `![Default]Value` column is
`z=log10(y/M)`. The script will do the inverse to generate the model
and pass the unit on to `.mod` files.

The Names of all tables must be unique.

### Compound

This table defines the compounds that are supposed to be modeled by
state variables and are subject to change by the reactions in the
systems. 

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| !Scale | log, log10, linear | and some variants of these|
| !InitialValue | a number | (per unit) in the above scale |
| !Unit | the unit of the above number | as it would be in linear scale |
| !SteadyState | `TRUE` | this compound should reach a steady state in at least one scenario and you want to know whether this happened |
|              |`FALSE` | it is not important whether or not this compound reaches steady state|
| !Assignment| `Name` or `ID`| this field will assign a pre-defined algebraic expression to the compound|

The conversion script will make a file called
`[…]SuggestedOutput.tsv`, it will have lines that can be used to check
whether a compound has reached steady state (or not), this is done for
each compound that has `!SteadyState` marked as `TRUE`. If that output
is close to `0`, then steady state was reached (it's the sum of all
fluxes for the compound in question).

Others columns are unused but may be informative to the user, or others.

#### Compound Assignments

In some cases, a compound's amount or concentration is not supposed to
be governed by reactions (kinetic laws, stoichiometry) but rather by a
fixed (time-dependent) value. In SBML this is called a boundary
condition. If the species is supposed to be constant, the field
`!IsConstant` can be set to `TRUE`; otherwise, you can assign the
value of an expression, listed in the `Expression` table to this
compound. The `!Assignment` field can contain the name of an
`Expression`. In SBML, a rule will be created in the `listOfRules`,
that rule will target the boundary condition species:

```xml
      <species    id="PKC_active" 
                name="PKC_active_value" 
         compartment="Comp1" 
initialConcentration="0" 
      substanceUnits="substance" 
   boundaryCondition="true"/>
<!-- ......... -->
<!-- and later -->
<!-- ......... -->
      <assignmentRule variable="PKC_active">
        <math xmlns="http://www.w3.org/1998/Math/MathML">
          <apply>
            <plus/>
            <ci> PKC_DAG_AA_p </ci>
            <ci> PKC_Ca_memb_p </ci>
            <ci> PKC_Ca_AA_p </ci>
            <ci> PKC_DAG_memb_p </ci>
            <ci> PKC_basal_p </ci>
            <ci> PKC_AA_p </ci>
          </apply>
        </math>
      </assignmentRule>

```

In the other formats, `vf` and `mod`, the relationship is much simpler:
```xml
<Expression Name="PKC_active_value" Description="defined expression Ex0" Formula="PKC_DAG_AA_p+PKC_Ca_memb_p+PKC_Ca_AA_p+PKC_DAG_memb_p+PKC_basal_p+PKC_AA_p"/>
<Expression Name="PKC_active" Description="defined expression S11" Formula="PKC_active_value"/>
```
which will lead to code such as:
```matlab
PKC_active_value = PKC_DAG_AA_p+PKC_Ca_memb_p+PKC_Ca_AA_p+PKC_DAG_memb_p+PKC_basal_p+PKC_AA_p;
PKC_active = PKC_active_value;
```

A blank cell, `NONE`, `FALSE`, or `NO` means that there is no
Assignment for this species/compound. Empty cells can be tricky if export or
import functions merge multiple delimiters (so `\t\t` is not
recognized as an empty cell).

A `TRUE` value in `!IsConstant` and meaningful assignments are
mutually exclusive and may lead to weird results.

### Parameters

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!Scale` | `log`, `log10`, `linear` | some aliases of these are possible (such as `base-10 logarithm`)|
| `!DefaultValue` | a number | in above scale, normalized to the unit of measurement, possibly subject to fitting/sampling |
| `!Std` | a number | standard deviation / uncertainty of this parameter |
| `!Min` and `!Max` | numbers | respectively, used if `!Std` is not present |

The columns `!Std` and `!Min`/`!Max` are only used in
sampling/optimization, the model conversion is unaffected by them, the
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

Since the kinetic law determines the reversibility of the reaction, the
column `!IsReversible` is not necessary, but if you determine the
kinetics based on the law of mass action it may be important for you
to have that column as a reminder (for when you are auto generating
the `!KineticLaw` column, which this script doesn't do).

### Input

The input parameters to the model that distinguish different
experiments. These quantities are known and can be influenced by the
people who are performing an experiment (or rather the real
counterparts of these quantities can be influenced). These play the
roles of (additional) parameters, but a different kind of parameter
than in the Parameter table. Experiments are supposed to have the same
parameters of the normal kind and different parameters of the input
kind.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!DefaultValue` | number | same as with normal parameters, but these are set to known values during an _experiment_ (these have to be known values) |
| `!ConservationLaw` | `TRUE`| the parameter is the result of an earlier run of the conversion script|
||`FALSE` | this parameter is unrelated to conservation laws |

The `!ConservationLaw` column will eventually be obsolete and the
numbers in that column will be determined by the _experiment-specific_
initial conditions. Currently there is only one _initial condition_ vector for
_all_ experiments.

At this moment, a conservation law parameter will be ignored whenever
the script runs (again), because the script creates those itself (every time it runs).

### Output

The outputs are _observable quantities_ of this system; what is and
isn't an output depends on what you can measure (or have knowledge
about). Outputs are usually converted to functions in the target
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
| `!ErrorName` | a string | indicates the column in data sheets that hold the measurement error of an observable |
| `!ErrorType` | not used | this is for the user |
| `!ProbDist`  | a string | the probability distribution of the noise model (currently unused); this is for humans to read |
| `!Formula`   | a math expression | the right hand side of the assignment |

Many data columns may share the same Error column. This is useful if
you have only a very rough estimate of the noise anyway, and the
outputs are in the same number range so using the same standard
deviation (etc.) for all data points seems good enough.

The data sheets are not used by this script, but they are used [here](https://github.com/a-kramer/mcmc_clib).

### Expression

These will be local variables of the model. Expressions are
assignments that are calculated repeatedly each time the ODEs right
hand side is called (before fluxes are calculated, fluxes are otherwise also a
type of _expression_).

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!Formula` | a math expression | right hand side of assignment|

Currently, there is no way to formally write something like _initial
assignments_, as those don't go into the model files, these
assignments have to happen before you call the model solver, so those
are up to you use of the model in the language you had in mind.

The only type of initial assignment that can be specified is the state
variables' initial conditions.

### Experiments

This table holds the mapping between input parameters and data
sheets. It determines the conditions under which a data set should be
replicated using the model.

| Column | Values  | Comment |
| -----: | :-----: | :------ |
| `!Type` | `Time␣Series` |  indicates that the data is a `t->output` mapping |
|       | `Dose␣Response` | data sheet is an input/output curve, i.e. `input->output` mapping |
| `>some_id` | a number | sets the input parameters for this experiment |
|`!Event`| a table name | the name of an event table that holds time instantaneous model state changes |

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
same sheet. Generally, it doesn't make sense to have the same ID for two 
rows, other than to override a general entry by a more specific entry.
