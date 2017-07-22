# SciRep2017
This folder contains the MATLAB functions and data used to generate all computational results presented in:

"Resolving the central metabolism of Arabidopsis guard cells"

Semidán Robaina-Estévez1*, Danilo de Menezes Daloso2,3*, Youjun Zhang2, Alisdair R. Fernie2, Zoran Nikoloski1**
1Systems Biology and Mathematical Modeling Group, 2Central Metabolism Group, Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany, 3Departamento de Bioquímica e Biologia Molecular, Universidade Federal do Ceará, Fortaleza, Brasil.

Scientific Reports, 2017


MainFunction contains the code to generate all results and calls the rest of the functions

GEM&DATA contains the Arabidopsis genome-scale model used in this study, the expression data (already preprocessed in R, see details in the publication) and two MATLAB structures: RESULTS and RESULTSCBC. These structures contain the solution generated with the G and M alternative optimal flux distributions sampled on this study when no additional constraints are imposed, and when a constraint on the ratio of carboxylation to oxygenation and reactions 11, 13 and 14 in the CBC are constraint to carry flux.The field "GCMlist" in each of these structures contains a list of all reactions in AraCOREred (split in forward and backward direction),  the mean flux values, the p-values of the Mann-Whitney tests and the maximum and minimum flux values calculated with RegrExFVA.

These functions depend on a working installation of Gurobi. Academic licenses can be requested from: http://www.gurobi.com"

Questions can be addressed to: robaina@mpimp-golm.mpg.de .In addition, the alternative optima space sample generated in this study can be requested under the same e-mail addresss.
