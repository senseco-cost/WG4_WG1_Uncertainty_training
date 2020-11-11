# WG4_Uncertainty_training

This is the code and data base for a course on uncertainty, calibration and traceability held on Romania at CETAL in fall 2019.


The structure and content is as follows:

- Schedule and Presentations: Includes the course schedule and the presentation that are explaining the tasks that are then implemented in the Matlab code.
- ExampleData: Data used by the provided Matlab code
- MatlabCode: Example Matlab code, covering the example and exercises of this training school

Matlab function details, data and relevant presentations follow below.



## Presentation:  Read the spectrometer data and determine statistics.pdf

** Data files:
| File name                                                   | Comment                                                                | 
| :---------------------------------------------------------- | :--------------------------------------------------------------------- | 
| 10000fL_15ms.xlsx
| 1000fL_15ms.xlsx
| 100fL_15ms.xlsx
| 5fL_15ms.xlsx
| Dark15ms.xlsx

** Matlab code: (note: the pathnames need adjusting to your local machine!)

ReadSpectrometerDataMeanStd.m


=======================
## Presentations: 

Intro to Monte Carlo and Application to RAD CAL.pdf
Radiometric Calibration Coefficient Determination.pdf


** Data files (see also last slide of presentation): 


L_Sphere.mat, STD_DN.mat, DN_L_CAL.mat, uL.mat		Input files for the code

u_rad_coeffs.mat		Output of Monte Carlo run: uncertainties of gain, offset and uncertainty due to gain and offset correlation



MC_Introduction.m		Code to produce plots shown in the intro to Monte Carlo

RAD_CAL_with_Linear_Fit_and_uncertainty_estimation_with_Monte_Carlo.m		Main script. Note: set the run_sim = true on line 408 to run MC (this can take very long!). Set to false once you have them calculated.
print_jpeg.m print_pdf.m		Functions to export figure to JPEG or PDF
progressbar.m gui_active.m		Functions for progress bar used to show progress during monte carlo run
get_realisations_gauss_dist.m	Function to create realisations



