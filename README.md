# Coupling Between NK Cell Cycle and Receptor Signaling

Paper Link:
https://scholar.google.co.in/citations?view_op=view_citation&hl=en&user=wxcnxHgAAAAJ&sortby=pubdate&citation_for_view=wxcnxHgAAAAJ:8k81kl-MbHgC

## Folder Descriptions 

### expt_data_analysis
This folder contains analyses related to Figure 2E for NKG2D-stimulated primary NK data (NKG2D_primary_NK_analysis.py).
Similar analyses were performed for NKG2D receptor-stimulated NKLs (NKG2D_NKL_analysis_G2_M_concatenated.py) and Ly49h-stimulated NKL cells (Ly49h_NKL_analysis_G2_M_concatenated.py).
The `data` subfolder contains the required data for the above analyses for primary NK cells and NKL cells (results are provided in the Supplementary Materials).
The `bootstrapping_shifted_geometric_mean` subfolder contains analyses of primary human NK cell datasets based on non-parametric tests, included in Supplementary Figure S12.

### ODE_vs_Gillespie_match
This folder contains the required files to generate Figures 3C-3D, demonstrating the match between the exact solution of ordinary differential equations and our in-silico Gillespie simulation for the linear protein model (Figure 3A).

### Gillespie_simulation
This folder contains the in-silico model for various cases discussed in Figures 4B-4C for the non-linear protein model (Figure 4A).

### Model_fitted_to_expt_data
In this folder, we fitted observed pS6 and pAkt datasets using our in-silico non-linear model (Figure 4D-4E) to determine whether the active or neutral scenario of receptor signaling better fits the data by estimating the phosphorylation rate.

### Gillespie_sensitivity_analysis_simulation
This folder contains sensitivity analyses examining how variations in the average abundance of a protein respond to a range of signaling parameters in our non-linear in-silico model (Figure 4A). Parameter ranges are discussed in the Supplementary Materials.