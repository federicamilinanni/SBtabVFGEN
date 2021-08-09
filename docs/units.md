# Units

SBML has support for units in all quantities. The units are specified using these 4 properties: `kind`, `scale`, `multiplier`, and `exponent` .
These unit attributes are interpreted like this: 

```
(multiplier * kind * 10^scale)^exponent
```
or, if you prefer: 
```
power(prod(multiplier,kind,power(10,scale)),exponent)
``` 

The `kind` can be any [SI](https://en.wikipedia.org/wiki/International_System_of_Units) base unit, e.g. `second`.
Other units can be derived using products of these base units, `liter/(nanomole millisecond)` is expressed as:

```
<unitDefinition id="liter_per_nanomole_millisecond">
 <listOfUnits>
  <unit kind="litre" exponent="1" scale="0" multiplier="1"/>
  <unit kind="mole" exponent="-1" scale="-9" multiplier="1"/>
  <unit kind="second" exponent="-1" scale="-3" multiplier="1"/>
 </listOfUnits>
</unitDefinition>

```

In _SBtab_ files, units are written in a human readable form (`!Unit`
column) and it is not always easy to interpret those units. The `R`
program in this repository attempts to read the units using a
regular expression with sub-groups. To make it somewhat doable, we have
additional rules on unit-strings:

1. Only SI base units are allowed for now
   - derived units such as `Newton` or `Hz` are not understood (liter is the only exception)
   - not all SI prefixes are understood, but the most common ones are (G,M,k,c,m,µ,n,p,f)
1. Only one slash is allowed (the slash has the lowest precedence)
   - `liter / mole second` is ok 
   - `1/((mole/liter) * second)` is not ok (because it has two slashes)
   - multiplication has higher precedence than division
1. All parentheses are ignored
   - `(mole/liter)^(-1)` is not interpreted correctly, because parentheses will be ignored
   - `liter/mole second` is the same as `liter/(mole second)`
   - `*` and `␣` are the same (blank space is interpreted as multiplication)
1. Powers cannot have spaces between base and exponent, the `^` is optional
   - `s^2` is ok
   - `s2` is ok and means the same thing
   - `kg m s^(-2)` is ok and the same as `kg m s^-2`
   - `kg m s-2` is also ok and means the same as `kg m s^-2`
   - `cm2` is ok and means square centimeters
   - `kg m s^( -2 )` is not ok (spaces)
1. The literal `1` is interpreted as: this quantity is dimensionless (in sbml this is actually called `dimensionless`)
   - a `1` in a unit will reappear in sbml, even if unnecessary
   - `1 m / 1 s` will have 4 entries in the sbml unit definition, with two unnecessary `dimensionless` units
   - `1/s` will be the same as `s^-1` in effect, but `1/s` will have an unnecessary `dimensionless` unit entry in the definition
   - no simplification of the unit is performed, so `meter/meter` is not simplified to `1`
1. long words can be used as well as abbreviations
   - `millisecond` is ok
   - `ms` is also ok
   - `msecond` (probably) also ok (but weird)
1. multipliers are always `1` (we don't have a system for
   multipliers). They are used to convert to and from non SI systems
   (imperial and so forth)
   - inches and feet are not multiples of powers of ten of SI base units.
   - `inches` are not parsed by the regular expression anyway (and we don't plan to ever support non SI units)
   - the only hope for _non SI_ units would be the _natural units_
1. The units that are not understood, but don't lead to a crash/error in the R code are interpreted as `1`
   - the user can then correct those definitions in the sbml file by hand (using a normal text editor [but please not notepad, treat yourself])

Because we use a quite simple regular expression to parse units, the base unit for mass is `g` (not `kg`, it's a bit of an exception in the SI world). Here is the regular expression, possibly not up to date:
```R
pat <- paste0("^(G|giga|M|mega|k|kilo|c|centi|m|milli|u|μ|micro|n|nano|p|pico|f|femto)?",
	          "(l|L|liter|litre|g|gram|mole?|s|second|m|meter|metre|K|kelvin|cd|candela|A|ampere)",
	          "\\^?([-+]?[0-9]+)?$")

```

   
## A list of perfectly fine units

Newton, Hertz and M are not parsed automatically, use the long form in column 2.

|meaning|suggested string| kind | scale | exponent |
|------:|:---------------|:----:|:-----:| :-------:|
|Newton| `kg m s-2` | `"gram"` | +3 | 1 |
||| `"metre"` | 1 | 1 |
||| `"second"` | 1 | -2 |
|nanomolarity (`nM`)| `nmol/l` | `"mole"` | -9 | 1 |
||| `"litre"` |  1 | -1 |
|kHz| `ms^-1`| `"second"` | -3 | -1 |

However, the units are not always interpreted right by importers (in other software).
Note that `kHz` is a bit unusual as Hertz is reciprocal to a base unit (second), so `kHz` is  `1000/s = 1/0.001 s = ms^-1`. 
