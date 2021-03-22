# VSAT-3D

VSAT-3D is part of the Valparaíso Stacking Analysis Tool (VSAT), it provide a series of tools for selecting, stacking, and analyzing 3D spectra. It is intended for stacking samples of spectra belonging to large extragalactic catalogs by selecting subsamples of galaxies defined by their available properties (_e.g. redshift, stellar mass, star formation rate_) being possible to generate diverse (_e.g. median, average, weighted average, histogram_) composite spectra. However, it is possible to also use VSAT on smaller datasets containing any type of astronomical object.

![Alt text](./Figures/Scheme1.jpg?raw=true "3D datacube Stacked spectra Scheme.")

## Content

1. .py:
   - Location of the input catalogue and spectral data. 
   - Parameters for selecting subsamples of galaxies according to their physical properties.
   - Stacking and Boootstrap parameters.
   - Location of the resulting products of the stacking analyses _e.g. tables, plots, pre-stacking processed spectra and stacked spectra_.

2. .py:
   - Math functions (e.g. cosmological constants, gaussian, lorentzian and voigt profiles for line emmision/absorption fitting) needed throughout the stacking analysis.

3. .py 
   - Functions to read, write and modify different tables. 

4. .py 
   - Funtions to access and modify (add, modify, delete) fits headers.

5. .py
   - Plot templates used throughout the stacking analysis. 



## Parameters
It is possible to ...

###### "Pre-1"
**pre_** 

###### "Pre-2"
**pre_smooth** (_bool, optional_) enables the spectral smoothing, **pre_smooth_shape** selects the smothing kernel (_i.e. gaussian,boxcar,mexican_) and **pre_smooth_size** sets the kernel size in pixel units.



![Alt text](./Figures/step.jpg?raw=true "Pre-processing of stacked spetra.")

## Line fitting
VSAT-3D uses lmfit for line fitting and ...

![Alt text](./Figures/FitSingle.jpg?raw=true "Pre-processing of stacked spetra.")


![Alt text](./Figures/FitMultiple.jpg?raw=true "Pre-processing of stacked spetra.")
## Example

The following example will stack a sample of 27 galaxies belonging to the Valpara\'iso ALMA/APEX Line Emission Survey(VALES). The sample of spectra can be downloaded from the 
[zenodo repository](). Then by simple running ```python Example.py``` will complete all the following steps below.

###### "Stacking"
The following snippet will stack galaxies from the COSMOS field. 

First . 

```python
os.clear()
```
This will generate the following tables saved in the tables directory (```~/Example/Stack_Results/COSMOS/TABLES/FRGRD```):
- P.csv

- 
Similarly ```Select_Subsamples``` can be used to define subsamples of objects restricted by any given object property. 

Next we stack the subsample. 

```python
os.clear()
```
This will generate the following fits files in the results directory (```~/Example/Stack_Results/COSMOS/STACKS/RESULTS/```):

```
 - P.fits
```

It will also generate ```*-c.fits```and ```*-smt.fits```files if ```pst_cnt``` and ```smt_spc_pst ``` are set to ```True```.

Auxiliary files generated during the stacking process are saved in ```~/Example/Stack_Results/COSMOS/STACKS/LAST-FITS/RESULTS```
and can be deleted. These files include a copy of the original spectra, continuum corrected, smoothed and interpolated spectra, continuum fit log files and .tex files with the spectra used in the stacking process. Notice that a same file can be used in differeent stackings and hence the -int.fits files include a number of the interpolated version.

###### "Stats"
Statistical values from the stacked galaxies can be obtained and saved in tables (```~/Example/Stack_Results/COSMOS/TABLES/STATS```) through:

```python
os.clear()
```

###### "Plots"
To plot the generated composite spectra (_e.g. average, median weighted average_)including post-processed versions (_continuum normalized and smoothed_) can bee generated with:

```python

```
This will generate and save a pdf plot file (```~/Example/Stack_Results/COSMOS/PLOTS/RESULTS/```) including the median, average, and weighted average in the upper panel and an histogram of the number of spectra combined per wavelength element in the lower panel.

![Alt text](./Figures/Stacked.jpg?raw=true "Stacked spectra computed COSMOS field.")

**plt_cnt_stk_spc** = ```True``` generates a plot of all the individual spectra used to generate the composite spectra (upper panel) with the composite specra in the bottom panel.

![Alt text](./Figures/Stacked-Contribution.jpg?raw=true "Stacked spectra COSMOS field.")

**plt_ind_spec** = ```True``` will generate individual plots (```~/Example/Stack_Results/COSMOS/PLOTS/IND//FRGRD/0-23/```) of all the spectra used to generate the compsite spectra. 

![Alt text](./Figures/Spec-Individual.jpg?raw=true "Stacked spectra COSMOS field.")

It will also generate a detailed plot of every step of the pre-processing procedure prior to the stacking process of all the combined spectra.
![Alt text](./Figures/Spec-Step.jpg?raw=true "Stacked spectra COSMOS field.")

###### "Line Fitting"

```python
os.clear()
```

All the line initial and fitted parameters (line center, amplitude, line width) will be saved in the corresponding composite fits header. For example for Lyman alpha the following headers will be created:

 - L05_CGLC=    1207.832434368174 / Lya1215.7 Ctr  1GF Crctlmfit

L05 corresponds to 

![Alt text](./Figures/LINE-FIT-COSMOS-avg-c-smt-G-Ind-Splt.jpg?raw=true "Stacked spectra COSMOS field.")	
```python
os.clear()
```
![Alt text](./Figures/LINE-FIT-COSMOS-avg-c-smt-G-Mlt-Splt.jpg?raw=true "Stacked spectra COSMOS field.")	

###### "Bootstrap"
To compute the Confidence Inteervals (CIs) of the composite spectra it is possible to bootstrap the spectra used in the stacking process. 

**bs_iteration_num** defines the number of bootstrap repetitionsm, **bs_percentage** defines the percentaje of elements to be resampled in each iteration (_default=1_), **bs_function**  defines the function used for the bootstrap.

It is possible to complete a broken BS process. 
**comp_BS_run**	(_bool_) enables the BS completion
**bst_brk_int** sets the reinizilation step at which the process broke, **bst_brk_lst** sets the BS number repetitions (_default=bs_iteration_num_)  #Last iteration step (default = bs_iteration_num) and **crts_BS_otb** (_bool_) creates the BS master table only for Completed BS Repetitions(_= bs_iteration_num_) but without master table. 


To generate the CIs, first we define:

```python
prefixF       = 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
tables        = [frg_dir_res + prefixF + '0-23.csv']
bs_function_s = [ '_med-c-smt']
```

Then we run the BS process:

```
os.clear()

```

This process will create into the ```~/BOOTSTRAP/STACKS/``` directory three diferent subdirectories: 
 - ```~/BOOTSTRAP/STACKS/LAST-FITS``` contains all the files used to generate stacked boootstrap in each repetition with a similar structure similar to  ```~/Example/Stack_Results/COSMOS/STACKS/RESULTS/``` directory and described above in the **Stacking** section.
 - ```~/BOOTSTRAP/STACKS/STATS-STR``` contains all the stacked boootstrap repetitions (_e.g. FILE_NAME-BS-1-stk-avg.fits, FILE_NAME-BS-2-stk-avg.fits, ... FILE_NAME-BS-N-stk-avg.fits_) stacked boootstrap repetitions.
 - ```~/BOOTSTRAP/STACKS/STATS-BST``` contains the bootsrap stacked spectra considering all the _N_ stacked boootstrap repetitions.

Finally, to plot the Stacked spectra including the CIs generated thorugh the BS process; first we define:

```
prefixF 	= 'P_Fg_' + CAT_PARENT + '_0-40-ss-zf_B-3-44-ss-zf_F-3-44-ss-sep_as-'
separation      = ['0-23'] 			#sep_as
bs_function_s   = ['med-c-smt']
```

and then

```
os.clear()
```

![Alt text](./Figures/P_Fg_COSMOS_BS_MST_100_med-c-smt-1150-1900.jpg?raw=true "Stacked spectra COSMOS field.")	

This plot will be saved in ```~/Example/Stack_Results/COSMOS/BOOTSTRAP/PLOTS/RESULTS/```.
**Notice that this plots is only useful to visualize the distribution of spectra generated through the bootstrap.** To properly generate CIs of the emission/absorption EW lines measured above, all the bootstrapped spectra should be fitted to obtain a distribution of EW. This can be easily done with:

```
os.clear()
```
If ```ivl_fts_hdr=True``` then the initial guess values will be read from the fits files headers, otherwise they will be set by the defualt line quantities defined in the Line_Dictionary.py file .
This will generate individual plots for each line profile fitted, located in ```~/Example/Stack_Results/COSMOS/BOOTSTRAP/PLOTS/INDIVIDUAL-BS/```.

As explained above all the fitted parameters (_i.e. line properties (amplitude, sigma, line center), fitted values_) are saved in the corresponding composite spectra fits headers.
## Dependencies
Currently VSAT works only with astropy 2.0 as it relies on pyraf continuum task for continuum normalization. However a new version will be released dropping this dependency.
 - [astropy](https://www.astropy.org)
 - [bottleneck](https://pypi.org/project/Bottleneck/)
 - [pandas](https://pandas.pydata.org)
 - [scipy](https://www.scipy.org)
 - [numpy](https://numpy.org)
 - [lmfit](https://lmfit.github.io/lmfit-py/)
 - [matplotlib](https://matplotlib.org)
 - [termcolor](https://pypi.org/project/termcolor/)
 - [progressbar](https://pypi.org/project/progressbar2/)
## License

BSD 3-Clause License

Copyright (c) 2021, VSAT-1D developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

