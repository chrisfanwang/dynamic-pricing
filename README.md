# Transfer Learning for Nonparametric Contextual Dynamic Pricing

This repository contains the implementation of the algorithms and simulations presented in the paper *Transfer Learning for Nonparametric Contextual Dynamic Pricing*. The code is designed to experiment with dynamic pricing strategies using transfer learning techniques to optimize revenue in markets with limited historical data.

## Overview
Dynamic pricing strategies adjust prices based on customer and market conditions to maximize revenue. This project focuses on leveraging historical data from related markets (source domain) to improve pricing decisions in a new market (target domain) under a covariate shift model. The Transfer Learning for Dynamic Pricing (TLDP) algorithm is introduced alongside existing pricing methods such as ABE and ExUCB.

## File Structure

### Core Algorithms
- **basic functions for TLDP.R** - Utility functions for TLDP, including computing ball partitions, selecting prices, and managing exploration.
- **dynamic pricing functions for auto loan data.R** - Implements TLDP and other pricing strategies specifically for the auto loan dataset.
- **dynamic pricing functions for Scenario_1.1.R** - Implements TLDP and benchmark methods for Scenario 1.1.
- **dynamic pricing functions for Scenario_1.2.R** - Implements TLDP and benchmark methods for Scenario 1.2.
- **dynamic pricing functions for Scenario_2.1.R** - Implements TLDP and benchmark methods for Scenario 2.1.
- **dynamic pricing functions for Scenario_2.2.R** - Implements TLDP and benchmark methods for Scenario 2.2.

### Scenario Experiments
- **Scenario_1.1.R** - Runs simulations varying parameter C_r.
- **Scenario_1.2.R** - Runs simulations varying parameter kappa.
- **Scenario_2.1.R** - Runs simulations varying parameter C_I.
- **Scenario_2.2.R** - Runs simulations varying the size of the source dataset n_P.

### Data Processing & Model Training
- **auto_loan.R** - Prepares and processes auto loan data, trains a random forest model, and applies TLDP to predict optimal prices.
- **Scenario_plots.R** - Aggregates and visualizes results from multiple scenarios using ggplot2.

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
1. Load required packages in R.
2. Run `auto_loan.R` to preprocess the dataset and train models.
3. Execute individual scenario scripts (`Scenario_*.R`) to run simulations and obtain results.
4. Run `Scenario_plots.R` to generate visualizations.

## Citation
If you use this code, please cite our paper:
> **Transfer Learning for Nonparametric Contextual Dynamic Pricing**

For any questions or issues, feel free to reach out!

