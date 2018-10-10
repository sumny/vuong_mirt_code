# vuong_mirt_code

This repository contains simulation code and results as well as other code to
replicate all findings reported in the paper [Model Selection of Nested and
Non-Nested Item Response Models using Vuong Tests]().

### Contents and getting started
* orig (directory containing the original simulation results)
* repl (directory containing possible rerurn simulation results)
* replication.R (R code for the simulations, text results and figures)

To inspect the original simulation results, change the working directory in R
to `vuong_mirt_code/orig`. To rerun the simulations, change it to
`vuong_mirt_code/repl` as described in `replication.R` and run the code in
`replication.R`.

### Packages and versions
* R version 3.4.3
* require(MASS) # version 7.3.48
* require(mirt) # version 1.26.3
* require(nonnest2) # version 0.5
* require(SimDesign) # version 1.7

### Authors
Lennart Schneider, R. Philip Chalmers, Rudolf Debelak, Edgar C. Merkle.

### License
This project is licensed under the terms of the GNU General Public License
v3.0.

