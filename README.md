gnominator
==============================
[![DOI](https://zenodo.org/badge/366433171.svg)](https://zenodo.org/badge/latestdoi/366433171)

Mixture model Naive Bayes Classifier for nominating representative microbiomes from two cohorts for gnotobiotic experiments

Naomi G. Wilson & Ariel Hernandez-Leyva - Conceptualization, Methodology, Software
S Joshua Swamidass, MD PhD - Conceptualization, Methodology
Andrew L. Kau, MD PhD - Conceptualization, Methodology, Supervision, Management, Funding Acquisition
@Washington University School of Medicine

See publication for more info: 
The gut microbiota of people with asthma influences lung inflammation in gnotobiotic mice.
Naomi G. Wilson, Ariel Hernandez-Leyva, Anne L. Rosen, Natalia Jaeger, Ryan T. McDonough, Jesus Santiago-Borges, Michael A. Lint, Thomas R. Rosen, Christopher P. Tomera, Leonard B. Bacharier, S. Joshua Swamidass, Andrew L. Kau
iScience 2023; doi: https://doi.org/10.1016/j.isci.2023.105991

## Installation Guide

built with R 3.6.3

Not working in R 4.x yet!

### Package dependencies

Users should install the following packages prior to installing `gnominator`, from an `R` terminal:

```
install.packages(c('ggplot2', 'caret', 'stats', 'pROC', 'RColorBrewer', 'reshape2'))
```

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('naomiwilson/gnominator', force=TRUE)
require(gnominator)
```

The package should take <30 seconds to install on a recommended computer. 

## System Requirements

### Hardware Requirements

The `gnominator` package requires only a standard computer with enough RAM to support the operations defined by a user. Authors used machines with these specs for development and calculating runtimes: 16 GB RAM, 6 2.2-Ghz cores.

### Software Requirements

#### OS Requirements

The developmental version of the package has been tested on the following systems:

macOSX 10.16

Windows 10

#### Runtimes
These runtimes are based on use of recommended machine with input of 95 microbiomes and 2125 unique ASVs total.

train_NBC: runtime 0.55 seconds

perform_NBC: runtime 1.176 seconds

