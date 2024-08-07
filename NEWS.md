# HydroPortailStats 1.1.0

* New Feature: estimation using historical flood data, supporting censored/interval data and systematic errors, as described in [Neppel et al., 2007](https://hal.science/hal-00509088). This is mostly a re-implementation of the Fortran `HBay` [executable](https://github.com/benRenard/BMSL/tree/main/cli/HBay), providing the following functions: `Hydro3_HBay()` to perform Bayesian estimation, `HBay_Plot()` to plot the result of the former function, and `Import_HBayConfig()` to import configuration files that were used by the `HBay` executable.
* New Feature: new distribution 'Triangle' available in distribution library
* New Feature: new function `GenerateWithinBounds()` to generate random values from a truncated distribution
* Bug fix: `Hydro3_Estimation()` was not handling p-to-T conversion properly when options$p2T<1.
* `HydroPortailStats` now has a logo

# HydroPortailStats 1.0.3

* Change impacting function interface: `GetQfromT(T,H3,options)` changed to `GetQfromT(RP,H3,options)` to avoid using T as a variable name
* Change impacting function interface: output of `GetTfromQ()` is now a list with names (RP,IC) rather than (T,IC), for consistency with the change above
* Updated HydroPortail URL in README and DESCRIPTION
* Miscellaneous changes for CRAN release with no impact on function interfaces
* Added a `NEWS.md` file to track changes to the package.

