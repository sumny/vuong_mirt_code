# Simulation code and results

This repository contains simulation code and results as well as other code to
replicate all findings reported in the paper [Model Selection of Nested and
Non-Nested Item Response Models using Vuong
Tests](https://doi.org/10.1080/00273171.2019.1664280)
([arXiv preprint](http://arxiv.org/abs/1810.04734)).

### Contents and getting started
* orig (directory containing the original simulation results)
* repl (directory containing possible rerun simulation results)
* replication.R (R code for the simulations, text results, tables and figures)

To inspect the original simulation results, change the working directory in R
to `vuong_mirt_code/orig`. To rerun the simulations, change it to
`vuong_mirt_code/repl`. Finally, run the code in `replication.R`.

### Packages and versions
* R version 3.5.1
* MASS version 7.3-51.1
* mirt version 1.29
* nonnest2 version 0.5-2
* SimDesign version 1.13

### Authors
Lennart Schneider, R. Philip Chalmers, Rudolf Debelak, Edgar C. Merkle.

### License
R code provided here is licensed under the terms of the GNU General Public
License v3.0.

