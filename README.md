# ContiBAIT

[![Build Status](https://travis-ci.org/oneillkza/ContiBAIT.svg?branch=master)](https://travis-ci.org/oneillkza/ContiBAIT)

Using strand inheritance data from multiple single cells from the organism whose genome is to be assembled, contiBAIT can cluster unbridged contigs together into putative chromosomes, and order the contigs within those chromosomes.

##Installation
Currently ContiBAIT is under review at Bioconductor, but once it is accepted, you can install it via:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite('ContiBAIT')
```

If you want the development version, you can get it from GitHub:

```{r}
# install.packages("devtools")
devtools::install_github("oneillkza/ContiBAIT")
```
