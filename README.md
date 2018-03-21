# Spatial Bayesian Latent Factor Model for Coordinate Based Meta-Analysis (SBLFRM-CBMA)

## Overview

This repository contains the Matlab code that runs the MCMC for the Spatial Bayesian Latent Factor Model described in Montagna et al. (2017).  Details of the analysis are outlined in Web Appendix B, and the data used for the neuroimaging application  is in Section 3 of the paper. 

## Reference

Montagna, S., Wager, T., Barrett, L. F., Johnson, T. D., & Nichols, T. E. (2017). Spatial Bayesian latent factor regression modeling of coordinate-based meta-analysis data. Biometrics. http://doi.org/10.1111/biom.12713

## Requirements

SPM software, available at http://www.fil.ion.ucl.ac.uk/spm/software/spm12/


## Input Data (3 files)

* `foci.txt` This file must contain the data of the neuroimaging application outlined in Section 3 of the paper. The data matrix has four columns: Study ID, X, Y, Z coordinates. The number of rows corresponds to the total number of foci. Coordinates are intended to be in the MNI coordinate space. Please convert coordinates measured in other spaces to MNI and make sure coordinates fall within a standard MNI brain mask. For example, to convert Talairach coordinates to MNI open Matlab and add the tal2mni function (see below) to the path, then run the following lines of code which will create the file `T88_new.txt`, containing the Talairach coordinates converted to MNI space.

 	  
```matlab
% File T88_old.txt contains the Talairach coordinates
T88 = textread('T88_old.txt'); 	
% addpath to spm (replace with your local SPM path)
addpath ~/Documents/spm
spm
T88_new = tal2mni(T88);
dlmwrite('T88_new.txt', T88_new);
```


* `covariates.txt` Thie file contains the study-specific covariates. The matrix has dimension N x r, where N is the number of studies and r is the number of covariates. In the paper, r = 5 and the covariates are modality (fMRI = 1 / PET = 0), the (standardised) number of subject scanned, inference method (fixed = 0 / random = 1 effects), type of stimulus (auditory = 0, visual = 1, recall = 2, imagery = 3, visual and auditory = 4, olfaction = 5), and p-values correction (corrected = 0 / uncorrected = 1)

* `studytype.txt` This file contains the study-types for the studies in the neuroimaging application in Section 3. The types are coded as 1 = Anger, 2 = Disgust, 3 = Fear, 4 = Happy, 5 = Sad

## Matlab code (4 files)

* `main_3D.m` This is the main code that runs the MCMC for the 3-dimensional spatial Bayesian latent factor model. Line-by-line explanation of each command is given in the file itself. It is important that the user sets the working directory to the source code and data location (modify line 4 as appropriate). Also, it is necessary to add the path to spm to use the nifti function at line 101. A 4 x 4 x 4 grid is used to approximate the integral at the HMC step (update of θ<sub>i</sub> in Web Appendix B). File mask_4mm.nii contains the brain mask used for the approximation. A different grid can be used, but will require modifications of the code at lines 53-57, 101, and 107. At line 100, the path to spm version 8 is added. If a different version of spm is used, please modify line 100 as appropriate. The code automatically saves the posterior samples for inference

* `HMC.m` This file includes the HMC function to update the basis function coefficients (Web Appendix B, update of θ<sub>i</sub>). This file should not be modified by the user

* `randraw.m` This file includes instrumental functions (it should not be modified by the user) 

* `tal2mni.m` This function converts coordinates to MNI brain best guess from Talairach coordinates (see above)

 
## R code (1 file)

* `Plots.R` This is the R code needed to reproduce some of the plots shown in the paper


## Brain Masks (2 files):

* `brainmask.nii` This is a 2 x 2 x 2 mm brain mask [91 x 109 x 91 voxels]

* `mask_4mm.nii` This is a 4 x 4 x 4 mm brain mask [45 x 54 x 45 voxels]
