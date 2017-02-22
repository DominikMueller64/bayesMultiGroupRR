# bayesMultiGroupRR

This package provides a [Stan](http://mc-stan.org/) model for fitting a bayesian hierarchical
multi-group regression model. For each group, a separate set of locus effects is estimated.
Across groups, locus effects are modeled to follow a multivariate normal distribution. 
## Installation

You can install bayesMultiGroupRR from github with:

```R
# install.packages("devtools")
devtools::install_github("DominikMueller64/bayesMultiGroupRR@master", build_vignettes = TRUE)
```

## Example

For an example of its usage, see

```R
library('bayesMultiGroupRR')
help(bayesMultiGroupRR)
```
