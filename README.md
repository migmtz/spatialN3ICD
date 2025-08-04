# spatialN3ICD

Source code for [Jagged-mediated lateral induction patterns Notch3
signaling within adult neural stem cell populations](link) [[1]](#1).

This Python code implements spatial analysis from point process theory to study Notch3 signalling between Neural Stem and Progenitor Cells (NSPCs) in the telencephalon of adult zebrafish.
Notch3 signaling is related to NSC stemness and quiescence states and the spatial correlation between cell nuclei and Notch3 activity levels are studied.

This study makes use of the marked K- and L-functions as described in [[2]](#2) along with statistical tests to detect deviations from a model of uncorrelatedness between spatial positions and mark values.

A detailed presentation of the methodology can be found in the Methods section of [[1]](#1) under the Spatial Statistics subsection.

## Data description

3 datasets are studied, corresponding to each folder.
Correspondances between folder names and lines or experiments in paper [[1]](#1) are as follows:
The Crispr folder corresponds to TgKI<sup>n3AG/+</sup>.
The BAC folder corresponds to TgBAC<sup>n3-GFP/+</sup>.
The Morpholino folder contains both control and treated group for fish injected with vivo morpholino against ligand Jagged1b.

## Group classification

In this paper, cells are often classified in groups corresponding to the value of their marks.
In general, by considering the quantiles of order q and 1-q, 3 groups are obtained classified in Lower, Medium and Upper groups.
Lower corresponds to cells with mark lower than quantile of order q.
Medium corresponds to cells with mark between quantiles of order q and 1-q.
Upper corresponds to cells with mark higher than quantile of order 1-q.

## Files description and image generation

* ```main_functions.py```: functions for computing both spatial and marked K-functions along with auxiliary functions.

Each folder may contain one of the following files:
* ```distribution_groups.py```: File to plot distribution of observed marks (m<sup>i</sup>) colored by groups. Generates Figure 3A.

* ```estimation_mK.py```: File to estimate the spatial and marked K-functions. Estimations are saved in saved_estimations using Pickle package.

* ```neighbourhoods.py```: File to plot the relative neighbourhoods of all observed cells in a radius r colored by groups. Generates Figures 3C, 6B, S4C.

* ```pairplot.py```: File to plot marks with respect to another reference mark in the observations. Regions are colored by groups and points are colored according to the intensity of the marks in the y-axis. Generates Figures S4A, S4B.

* ```plot_mk2_tests.py```: File to compute and plot the deviation test from an uncorrelated model using the estimated marked K-functions. Also used to compute the corresponding p-values. Generates Figures 3D, 6C, S4D.

Each folder contains two folders: ```images```containing the images found in the paper, ```saved_estimations``` for the saved estimated K-functions.

## References

<a id="1">[1]</a>
Ortica, S., Martinez Herrera, M., Degroux, L., Rochette, B., Dray, N., Bally-Cuif, L., Jagged-mediated lateral induction patterns Notch3
signaling within adult neural stem cell populations (2025). bioRxiv: [⟨10.1101/2025.07.29.667421⟩]([link](https://www.biorxiv.org/content/early/2025/08/01/2025.07.29.667421)).

<a id="2">[2]</a>
Penttinen, A., Stoyan, D., and Henttonen, H.M. (1992). Marked Point Processes in Forest Statistics.
Forest Science 38, 806–824. [doi](https://doi.org/10.1093/forestscience/38.4.806).
