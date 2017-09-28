# Pulse stream models in ultrafast ultrasound imaging
[Adrien Besson](https://adribesson.github.io/), [Rafael E. Carrillo](https://scholar.google.com/citations?user=-KhnCwMAAAAJ&hl=fr), [Dimitris Perdios](https://people.epfl.ch/dimitris.perdios?lang=fr), [Marcel Arditi](https://scholar.google.ch/citations?user=4w3-BxEAAAAJ&hl=fr), [Yves Wiaux](https://researchportal.hw.ac.uk/en/persons/yves-wiaux) and [Jean-Philippe Thiran](http://lts5www.epfl.ch/thiran.html).
submitted to [IEEE ICASSP 2018](https://2018.ieeeicassp.org/default.asp) 

## Introduction
This repository contains all the code to reproduce the results of the paper: Pulse stream models in ultrafast ultrasound imaging.

## Abstract:
This paper considers the problem of reconstructing ultrafast ultrasound (US) element raw-data from random projections. It presents a new signal model, coined as multichannel ultrasound pulse-stream model, which exploits the pulse streams models of US signals and accounts for the intersensors dependencies. We propose a sampling theorem and a reconstruction algorithm, based on `1-minimization, for signals belonging to such a model. We show the benefits of the proposed approach through numerical simulations on 1D-signals and on in vivo carotid images.

## Citation (BibTex):
If you are using this code, please cite the following paper. 

## System Requirements:
This software has been tested on Linux 64-bit system (Ubuntu 16-04, Mint distribution) and on Windows 10 system.

### Prerequisites
1. MATLAB (Tested on R2017a)

## Running USStream:
1. Download and unzip USStream-master.zip
2. Open MATLAB and navigate in the folder USStream-master
3. If you want to reproduce the Figures of the paper, run the script'Script_reproduce_figure_SPL.m'
4. If you want to regenerate:
	* The results of the noiseless experiments (used to generate Figure 2), run the script 'Script_reproduce_noiseless_experiment.m'
	* The results of the noisy experiments (used to generate Figure 3), run the script 'Script_reproduce_noisy_experiment.m'

## Contact:
Adrien Besson, (adrien.besson@epfl.ch)

