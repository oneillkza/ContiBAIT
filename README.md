# ContiBAIT

[![Build Status](https://travis-ci.org/oneillkza/ContiBAIT.svg?branch=master)](https://travis-ci.org/oneillkza/ContiBAIT)

Using strand inheritance data from multiple single cells from the organism whose genome is to be assembled, contiBAIT can cluster unbridged contigs together into putative chromosomes, and order the contigs within those chromosomes.

## Installation

ContiBAIT is best installed via Bioconductor:

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('ContiBAIT')
```

If you want the development version, you can get it from GitHub:

```{r}
install.packages("devtools")
devtools::install_github("oneillkza/ContiBAIT")
```
