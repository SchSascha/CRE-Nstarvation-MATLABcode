# CRE-Nstarvation-MATLABcode
This code goes alongside the publication de Lomana et al. 2015 "Transcriptional program for nitrogen starvation-induced lipid accumulation in Chlamydomonas reinhardtii". The code should enable FBA based simulation of TAG and biomass differences upon N starvation.

A good start is calling/uncommenting a line of runAll for different experiments.
This script forwards calls to further functions and scripts to analyse, identify, and visualise biomass and TAG yield differences upon nitrogen starvation.

The individual models are created on the fly, depending on which PRISM reaction is used as light source and which biomass function is used.

PRE:
  running COBRA toolbox in MATLAB
