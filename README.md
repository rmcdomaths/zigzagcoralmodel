# Zig-zag persistence of a stochastic, spatial model for coal reefs

## Overview

This repository contains the code used in https://arxiv.org/abs/2209.08974. Data simulated for the paper may be simulated directly using the agent-based model https://github.com/rneuhausler/coralModel-TDA-study or downloaded from https://figshare.com/articles/dataset/zigzagcoralmodel_simulated_data/20409063. To simulate the data, clone the coralModel repo and run the generate_figX.sh files (X = 1,2,3,4,5). Note that simulation will take a few hours on standard desktops.  If downloading the simulated data from figshare instead, copy the folder named 'output' into this repo.

## TDA

Computation of Persistent Homology and Zig-zag persistence relies on the BATS package: https://github.com/CompTop/BATS.py. Once BATS has been installed, most functions for computing persistence (of images, using cubical complexes) may be found in TDAtools.py. Much of TDAtools.py is specialised for outputs of the coral model (see link above) but its functions may be adapted for the computation of zig-zag of other datasets. The files generate_bars_figX.sh files (X = 2,4,5) need to be run to generate zigzag bars and landscapes.

## Figures

The figures from the paper may be generated using the various .ipynb files. Run these with jupyter notebook. Note that the data folder ('output') must be present for this to work.
