# Pulse-stream models in time-of-flight imaging
[Adrien Besson](https://adribesson.github.io/), [Dimitris Perdios](https://people.epfl.ch/dimitris.perdios?lang=fr), [Marcel Arditi](https://scholar.google.ch/citations?user=4w3-BxEAAAAJ&hl=fr), [Yves Wiaux](https://researchportal.hw.ac.uk/en/persons/yves-wiaux) and [Jean-Philippe Thiran](http://lts5www.epfl.ch/thiran.html).
submitted to [IEEE ICASSP 2018](https://2018.ieeeicassp.org/default.asp) 

## Introduction
This repository contains all the code to reproduce the results of the paper: Pulse-stream models in time-of-flight imaging.

## Abstract:
This paper considers the problem of reconstructing raw signals from random projections in the context of time-of-flight imaging with an array of sensors. It presents a new signal model, coined as multichannel pulse-stream model, which exploits pulse-stream models and accounts for additional structure induced by inter-sensor dependencies. We propose a sampling theorem and a reconstruction algorithm, based on l1-minimization, for signals belonging to such a model. We show the benefits of the proposed approach by means of numerical simulations and on a real non-destructive-evaluation application.

## Citation (BibTex):
If you are using this code, please cite the following paper. 

## System Requirements:
This software has been tested on Linux 64-bit system (Ubuntu 16-04, Mint distribution) and on Windows 10 system.

### Prerequisites
1. MATLAB (Tested on R2017a)

## Running USStream:
1. Download and unzip USStream-master.zip
2. Open MATLAB and navigate in the folder USStream-master
3. If you want to reproduce the Figures of the paper, run the script 'Script_reproduce_figure_SPL.m'
4. If you want to regenerate:
	* The results of the noiseless experiments (used to generate Figure 2), run the script 'Script_reproduce_noiseless_experiment.m'
	* The results of the noisy experiments (used to generate Figure 3), run the script 'Script_reproduce_noisy_experiment.m'
	* The results of the noisy experiment (used to generate Figure 4), run the script 'Script_reproduce_invivo_experiment.m'
5. If you are interested in testing the least-square solution of the problem, explained in Section 2.4.1 of the paper, for the noiseless experiment, you can run the script 'Script_least_squares_noiseless_experiment.m'

## Contact:
Adrien Besson, (adrien.besson@epfl.ch)

