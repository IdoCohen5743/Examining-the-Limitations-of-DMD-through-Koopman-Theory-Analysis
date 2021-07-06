# Examining the Limitations of DMD through Koopman Theory Analysis

Code repository for our version of ["Examining the Limitations of DMD through Koopman Theory Analysis"](Essay_full_link) by [Ido Cohen](https://idoc.webgr.technion.ac.il/),  [Guy Gilboa](https://guygilboa.net.technion.ac.il/), at Technion - Israel Institute of Technology.

This work binds the existence of \ac{KEF}, the geometric of the dynamics, and the validity of \ac{DMD} to one coherent theory. Viewing the dynamic as a curve in the state-space allows us to formulate an existence condition of \acp{KEF} and their multiplicities. These conditions lay the foundations for system reconstruction, global controllability, and observability for nonlinear dynamics. \\
\ac{DMD} can be interpreted as a finite dimension approximation of \ac{KMD}. However, this method is limited to the case when \acp{KEF} are linear combinations of the observations. We examine the limitations of \ac{DMD} through the analysis of Koopman theory. We propose a new mode decomposition technique based on the typical time profile of the dynamics. An overcomplete dictionary of decay profiles is used to sparsely represent the dynamic. This analysis is also valid in the spatial continuous setting of Koopman theory, which is based on variational calculus.\\
We demonstrate applications of this analysis, such as finding \acp{KEF} and their multiplicities, calculating \ac{KMD}, dynamics reconstruction, global linearization, and controllability.


<p align="center">
	<img src="https://i.imgur.com/2LTDu9w.png" | height=250>
</p>



## Requirements
This code has been implemented using Matlab-2018b

## Data & directory tree

```
code_submission 	# code and data together
├── fast_experiment1D_final.m  		 # Full experiment script for toy example (Three pulses)
├── fast_experiment2D_zebra_final.m  # Full experiment script for harder example (Line from zebra image)
├── proj_tvl2.m          			 # Helper function for ss_freq_tv_evolve
├── ss_freq_tv_evolve.m   			 # Mathematical Engine of TV - evolve TV flow
├── TV_data1D.mat         		     # Data for toy example (three pulse) experiment
├── TV_ZebraLine.mat      			 # Data for zebra experiment
├── zebra_media_gmu.jpg   			 # Image from which the input signal for zebra experiment was taken
├── Additional Data       			 # More Data can be added here
```

## Run algorithm
```bash
# Run toy example experiment
fast_experiment1D_final
```

```bash
# Run zebra row experiment
fast_experiment2D_zebra_final
```
Full computational mode:
```bash
# In order to run either code without relying on pre-saved data in .mat files, please change the loadData flag in line 3 to 0
loadData = 0;
```

## Runtime
Running in full computational mode (without loading data) might change from machine to machine. Results presented for running experiment on a 8-Gen. Core i7 laptop, with 16BG of RAM.

- A full run of the basic toy example experiment takes ~8.5 Minutes
- A full run of the advanced experiment with the row from the zebra image takes ~2 Hours

The above mentioned long times are due to usage of standard method for computing subgradients iteratively, which is the method we compare to.

- All other run modes take less than 1 minute


## Citation
If you find our work useful, please cite our paper:
```bash
@InProceedings{10.1007/978-3-030-75549-2_5,
author="Cohen, Ido
and Berkov, Tom
and Gilboa, Guy",
editor="Elmoataz, Abderrahim
and Fadili, Jalal
and Qu{\'e}au, Yvain
and Rabin, Julien
and Simon, Lo{\"i}c",
title="Total-Variation Mode Decomposition",
booktitle="Scale Space and Variational Methods in Computer Vision",
year="2021",
publisher="Springer International Publishing",
address="Cham",
pages="52--64",
isbn="978-3-030-75549-2"
}


```

Feel free to place issues here or contact me via the e-mail in my personal page.
 
