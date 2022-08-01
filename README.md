# Zig-zag persistence quantifies the resilience of spatial patterns in a coral reef model

## Overview

This repository contains the code used in (paper link). Data simulated for the paper may be simulated directly using the agent-based model https://github.com/rneuhausler/coralModel or downloaded from https://figshare.com/articles/dataset/zigzagcoralmodel_simulated_data/20409063. To simulate the data, clone the coralModel repo and copy the .sh files into the main folder and run. Note that simulation will take a few hours on standard desktops.  If downloading the simulated data instead, copy the folder named 'output' into this repo.

## TDA

Computation of Persistent Homology and Zig-zag persistence relies on the BATS package: https://github.com/CompTop/BATS.py. Once BATS has been installed, most functions for computing persistence (of images, using cubical complexes) may be found in TDAtools.py. Much of TDAtools.py is specialised for outputs of the coral model (see link above) but its functions may be adapted for the computation of zig-zag of other datasets.

## Figures

The figures from the paper may be generated using the various .ipynb files. Run these with jupyter notebook. Note that the data folder ('output') must be present for this to work.
