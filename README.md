# SOM

This is an attempt to actually keep track of various updates and changes to the Statistical Oxidation Model, and to make it more available for use by others.

Questions should be directed to Chris Cappa (CDC; cdcappa@ucdavis.edu)
  - please report any bugs found; no guarantee is made of bug-free code
  - if you want to contribute updated code, please contact CDC first
  - if you use the SOM in your work, please let CDC know

Consists of five files:
1. SOM OH reactions - vX.Y.Z.pxp --> the main SOM program. The following four .ipf files should be embedded already, but are contributed for posterity and tracking
2. SOM procs.ipf --> various procedures associated with running the single-component SOM
3. SOM-MC procs.ipf --> various procedures associated with running the multi-component SOM
4. SOM fitting.ipf --> various procedures assocaited with fitting data with SOM
5. SOM graphs.ipf --> various procedures associated with recreating/creating graphs/panels/tables

To run SOM and have things compile properly, you will need the following additional .ipf files (available at https://github.com/chriscappa/Igor-Functions) 

1. GlobalUtils_cdc.ipf
2. 2015GeneralMacros.ipf
3. SizeDistributionProcessing.ipf
