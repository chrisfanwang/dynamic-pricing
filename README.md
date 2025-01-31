# Transfer Learning for Nonparametric Contextual Dynamic Pricing

This repository contains the algorithms and numerical experiments presented in the paper *Transfer Learning for Nonparametric Contextual Dynamic Pricing*. 

## Overview
This project focuses on leveraging source data to improve pricing decisions for target data under a covariate shift model. The Transfer Learning for Dynamic Pricing (TLDP) algorithm is introduced and compared with existing pricing methods including ABE and ExUCB.

## File Structure


- **basic functions for TLDP.R** - Utility functions for TLDP.
- **dynamic pricing functions for auto loan data.R** - Implements TLDP, ABE and ExUCB for the auto loan dataset.
- **dynamic pricing functions for Scenario_1.1.R** - Implements TLDP, ABE and ExUCB for Scenario 1.1.
- **dynamic pricing functions for Scenario_1.2.R** - Implements TLDP, ABE and ExUCB for Scenario 1.2.
- **dynamic pricing functions for Scenario_2.1.R** - Implements TLDP, ABE and ExUCB for Scenario 2.1.
- **dynamic pricing functions for Scenario_2.2.R** - Implements TLDP, ABE and ExUCB for Scenario 2.2.


- **Scenario_1.1.R** - Runs simulations for Scenario 1.1.
- **Scenario_1.2.R** - Runs simulations for Scenario 1.2.
- **Scenario_2.1.R** - Runs simulations for Scenario 2.1.
- **Scenario_2.2.R** - Runs simulations for Scenario 2.2.

- **auto_loan.R** -  Processes auto loan data and applies dynamic pricing algorithms.
- **Scenario_plots.R** - Aggregates and visualizes results using ggplot2.

## Dependencies
This project requires R and the following R packages:
- `dplyr`
- `ggplot2`
- `randomForest`
- `egg`
- `latex2exp`
- `cowplot`
- `patchwork`

## Usage
 Load required packages in R.
 Run `auto_loan.R` to preprocess the auto loan dataset. 
 Execute individual scenario scripts (`Scenario_*.R`) to run simulations and obtain results.
 Run `Scenario_plots.R` to generate visualizations.



