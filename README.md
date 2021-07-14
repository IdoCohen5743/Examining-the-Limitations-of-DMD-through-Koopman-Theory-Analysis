# Examining the Limitations of DMD through Koopman Theory Analysis

Code repository for our version of ["Examining the Limitations of DMD through Koopman Theory Analysis"](Essay_full_link) by [Ido Cohen](https://idoc.webgr.technion.ac.il/),  [Guy Gilboa](https://guygilboa.net.technion.ac.il/), at Technion - Israel Institute of Technology.

This work binds the existence of \ac{KEF}, the geometric of the dynamics, and the validity of \ac{DMD} to one coherent theory. Viewing the dynamic as a curve in the state-space allows us to formulate an existence condition of \acp{KEF} and their multiplicities. These conditions lay the foundations for system reconstruction, global controllability, and observability for nonlinear dynamics. \\
\ac{DMD} can be interpreted as a finite dimension approximation of \ac{KMD}. However, this method is limited to the case when \acp{KEF} are linear combinations of the observations. We examine the limitations of \ac{DMD} through the analysis of Koopman theory. We propose a new mode decomposition technique based on the typical time profile of the dynamics. An overcomplete dictionary of decay profiles is used to sparsely represent the dynamic. This analysis is also valid in the full continuous setting of Koopman theory, which is based on variational calculus.\\
We demonstrate applications of this analysis, such as finding \acp{KEF} and their multiplicities, calculating \ac{KMD}, dynamics reconstruction, global linearization, and controllability.


<p align="center">
	<img src="https://i.imgur.com/2LTDu9w.png" | height=250>
</p>



## Requirements
This code has been implemented using Matlab-2018b

## Data & directory tree

```
code_submission 	# code and data together
├── LinearDecay_EigenFunctional.m  		 # Full experiment script for toy example (Two functional)
├── mexLasso.mex  						 # Lasso implementation of from [1]
├── write_pdf_New_Image.m          		 # PDF creator
[1] M. Elad,Sparse and redundant representations:  from theory to applications in signal and image process-ing, Springer Science & Business Media, 2010
```

## Run algorithm
```bash
# Run toy example experiment
LinearDecay_EigenFunctional.m
```

## Runtime
Not enough to make a coffee


## Citation
If you find our work useful, please cite our paper:
```bash



```

Feel free to place issues here or contact us via the e-mail in my personal page.
 
