#faoutlier

Factor analysis outlier and influential case detection in R. Provides 
tools useful for detecting and summarize influential cases that
can affect exploratory and confirmatory factor analysis models. 

## Note
This package has two branches, one for CRAN (master) and one that utilizes 
[OpenMx](http://openmx.psyc.virginia.edu/ "OpenMx homepage") (omxversion). The OpenMx version can be installed in `R` by following these steps:

1. If not installed, obtain the `devtools` package with `install.packages('devtools')`
2. Install the github maintained `faoutlier` package using `devools::install_github('faoutlier', username='philchalmers', branch='omxversion')`. This will overwrite any previous version that may have been installed (and can not be automatically updated).
3. Now load the new faoutlier normally with `library(faoutlier)`

