gecko
==============================

Mixture model Naive Bayes Classifier for nominating representative microbiomes from two cohorts for gnotobiotic experiments

## Installation Guide

built with R 3.6.3

Not working in R 4.x yet!

### Package dependencies

Users should install the following packages prior to installing `gecko`, from an `R` terminal:

```
install.packages(c('ggplot2', 'caret', 'stats', 'pROC', 'RColorBrewer', 'reshape2'))
```

### Package Installation

From an `R` session, type:

```
require(devtools)
install_github('naomiwilson/gecko', force=TRUE)
require(gecko)
```

The package should take <30 seconds to install on a recommended computer. 

## System Requirements

### Hardware Requirements

The `gecko` package requires only a standard computer with enough RAM to support the operations defined by a user. Authors used machines with these specs for development and calculating runtimes: 16 GB RAM, 6 2.2-Ghz cores.

### Software Requirements

#### OS Requirements

The developmental version of the package has been tested on the following systems:

macOSX 10.16
Windows
