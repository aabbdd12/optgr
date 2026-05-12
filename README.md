# optgr: Optimal Group Targeting for Poverty Reduction

This repository hosts the official Stata package for **Optimal Group Targeting for Poverty Reduction (OPTGR)**, developed by Araar Abdelkrim.

## Description
The `optgr` command (via the `ogtpr` syntax) implements the sequential data-graph algorithm of Araar and Tiberti (2017) for optimal anti-poverty group targeting under a fixed per-capita budget. 

### Key Features:
* **Algorithm exactitude:** Characterises the small-transfer regime where classical water-filling is exact, and the large-transfer regime where the data-graph algorithm is required.
* **Vectorised performance:** Proposes a vectorised implementation achieving a 2.5 to 4.8x speedup.
* **Universal FGT Support:** Solves the allocation problem for all additive Foster-Greer-Thorbecke (FGT) poverty indices, including an adaptive partitioning extension to resolve the double-inflection problem for $\alpha=0$.

## Installation
You can install the latest version directly from this GitHub repository by typing the following command in Stata:

```stata
net install optgr, from("[https://raw.githubusercontent.com/Araar-Abdelkrim/optgr/main](https://raw.githubusercontent.com/Araar-Abdelkrim/optgr/main)") replace
```

## Datasets Included
The package includes the replication dataset **bkf98I.dta** (Enquête Prioritaire II of Burkina Faso, 1998) used in the methodological paper. Once installed, you can load it directly in Stata.

## Usage
To run the included example script:
```stata
do example.do
```

Basic usage syntax:
```stata
use bkf98I.dta, clear
ogtpr exppc, hgroup(gse) hsize(size) alpha(0) pline(80000) trans(4000) ered(1)
```

## Documentation
* For syntax details, run `help ogtpr` in Stata after installation.
* For the full theoretical background, please refer to the included working paper: `beyond_datagraph_FINAL.pdf`.