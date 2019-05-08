GAIL
====

Geo-Allocation of Irregular Units

This is an [R](https://cran.r-project.org/) package for allocating data from irregular spatial units onto a regular set of spatial units.

The motivating use-case for this package wasto allocate cases of some disease (e.g., influenza) from PO Box ZIP codes to nearby Standard ZIP codes. Often data from PO Box ZIP codes are discarded, or all of the data are aggregated to a larger spatial unit such as county. To preserve the cases associated with PO Box ZIP codes at the more fine spatial resolution (Standard ZIP code), we stochastically allocate cases from PO Box ZIP codes to 

Notes about the package:

- It is very new. There may be incomplete entries or placeholder text in some areas.
- The **sf** package is used to perform spatial operations. Since **sf** is also relatively new, **GAIL** is subject to change based on updates to **sf**.
- There are currently there are no data provided with the package to use for demonstration. Eventually a dataset will be provided for this purpose.

## Installation

**GAIL** is not yet on CRAN (and should mature for some time before submission is even considered).

To install, you can use:

`devtools::install_github( "jelsema/GAIL" )`

A Windows version may be added here in the future.

## Usage

Some basic examples will be added here in the future. Eventually a vignette will be provided with the package.


## License

This package is provided under the GPL-3 license.
