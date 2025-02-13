# MFPCA with categorical trajectories in R

This repository contains **R scripts** for performing **Multivariate Functional Principal Component Analysis (MFPCA)** on temporal categorical data. 
The goal is to apply MFPCA techniques to analyze and reduce the dimensionality of categorical datasets within a functional framework.

## Repository Contents

- **R Scripts**: 'TDS_MFPCA_2502.r and TCATA_MFPCA_250211.r. They allow you to reproduce the paper results.
- **Utility Functions**: Some functions used in the R scripts to preprocess the data that should be sourced at the beginning of the script.


## Installation

To use the scripts in this repository, follow these steps:

### Prerequisites

1. Make sure you have **R** installed. You will also need to install the following packages:

- `cfda`
- `ggplot2`
- `openxlsx`
- `gridExtra`
- `funData`
- `MFPCA`

You can install the required packages by running the following command in R:

```r
install.packages(c("cfda","ggplot2","openxlsx","gridExtra","funData","MFPCA"))
```

2. Clone this repository into your working environment:
```bash
git clone https://github.com/Chemosens/ExternalCode/new/main/MFPCAWithCategoricalTrajectories


 
