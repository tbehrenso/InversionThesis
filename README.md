# InversionThesis

Just a quick overview of the pipeline and which files it uses because its not that obvious and there are many redundant files.

## Scenario Naming
The naming scheme of the different scenarios is given below, and is used for labelling the SLiM scripts and any output files.

- ***neutral_2pop***: Completely neutral
- ***locallyAdapted_2pop***: No inversion, but two locally adapted alleles are intrdouced
- ***adaptiveInversion_2pop***: Inversion which is itself adaptive, but no locally adapted alleles
- ***inversionLAA_2pop***: Introduces an inversion containing two locally adapted alleles

## Pipeline Summary

1. **SLim Simulations**: Runs on the cluster controlled by [Inversion_Rep_Burnin.sh](Inversion_Rep_Burnin.sh), which sets most of the parameters and loops over a given SLiM script to create multiple MS outputs at different timepoints
    - Running Locally: Can also run each of the SLiM scripts locally (ex. [inversionLAA_2pop.slim](slim_scripts/inversionLAA_2pop.slim)), but need to un-comment lines that set the parameters within the script, and comment out any lines that save MS output to a specific directory.
2. **Calculations and Plotting**: On the cluster, is controlled by [R_Plotting.sh](R_Plotting.sh), which uses [Inversion_Rep_General.R](R_Scripts/Inversion_Rep_General.R) to do all the calculations and create the various plots.
    - Running Locally: [Inversion_Rep_General.R](R_Scripts/Inversion_Rep_General.R) should be able to run locally, as long as the working directory and PATH on lines 8 and 34 are updated. The script "detects" whether its running on the cluster based on the operating system (so Linux is assumed to be on the cluster). Not sure if MacOS would be interpreted as Linux or not by this function.
    - Some R scripts are not included in the main script (ex. the correlation heatmap is its own script: [justCorr.R](R_Scripts/justCorr.R)). For these, I use the [justOne_Plotting.sh](justOne_Plotting.sh) script, and adapt the settings/values accordingly.
    
## The R Script - Overview
My script is somewhat messy at the moment, but this is a quick overview of how its organized under the different sections
- **Prep**: The first few lines serve to load packages and determine some basic parameters
- **Parameters**: Set constants, and extract relevant parameters (depending on if cluster or local)
- **Functions**: A collection of functions. Some of which not used in the current code, but can give other information such as average inversion frequencies
- **Data Extraction**: Extracts *most* of the relevant data in a single loop through all the output files in the given directory. For the most part, this will have already put the data into sliding windows so it can be placed in a single dataframe per summary statistic.
    - The idea here was to only need to access each file once cuz this seems to be somewhat slow. But I end up having multiple other places in the code that loop over these same files. In some cases this is because I was lazy, and other times because some calculations require a different type of loop that acccesses both populations of a simulation at once. This does create the problem that code cannot be run for JUST one summary statistic, so I usually end up copying the relevant sections of the loop to a separate script for this purpose.
- **Diversity**: Calculations of both Hexp and Nucleotide Diversity. This includes calculation for nucleotide diversity that separates haplotypes, which is in its own loop.
- **Correlation**: Averages across the 3D matrix (# of windows x # of windows x # of replicates) created from the data extraction step for each population and plots as heatmap
- **Differentiation**: Calculations for Fst between populations and between haplotypes, each in their own separate loop (fully independent of data extraction loop). The calculations for Fst between haplotypes are skipped when the inversion shouldn't be present.
- **Finalization**: Outputting plots into graphics device or as a file, depending on where script is being run.
    - Note: Currently, any plots made using ggarrange() do not show up in R-Studio when sourcing the script, but do if they are run manually.
