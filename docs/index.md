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

## SBtab

This format is human readable and writable. We use the standard as we understand it, but use our own code to process the files, using a [subset of named columns](./sbtab.md). These files can hold a great variety of content. 
