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

First, we stack the subsample. 

```python
for element in itlpd(channel_width,sbsmn,sbsms):

	nchan         = (2*subcube_width /element[0])+1
	slice_nmb     = int(np.ceil(((nchan-1)/2)))#-1 #19 
	#Results
	stk_hme_dir  = home + 'Stack_Results-'+ line +'-3D/'
	img_dir_res  = stk_hme_dir + 'IMAGES/'  #+ str(element[0]) +'/'
	stp_dir_res  = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res  = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_res  = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stk_dir_res  = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'
	DIR_CAT_IPT  = [cats_dir]
	DIR_SPC_IPT  = [img_dir]
	DIR_RES      = [stk_hme_dir,stp_dir_res,tbl_dir_res,plt_dir_res,stk_dir_res]
	Check_directories(cat_ipt_tbl,cat_parent)

	cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

	print colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan')
	print

	Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
	fits          = Cat_Ipt_Tbl[2]
	delta_nu      = Cat_Ipt_Tbl[4]
	z             = Cat_Ipt_Tbl[8]  #[5]
	Lfir          = Cat_Ipt_Tbl[11] #[6]
	nu            = Cat_Ipt_Tbl[13] #[14]
	vel           = Cat_Ipt_Tbl[14] #[15]
	num_obj       = len(Cat_Ipt_Tbl[0])

	cubetoread = [(img_dir_res + image_fits+ '.'+str(element[0]) +'kms.fits') for image_fits in (fits)]
	print colored('Reading files as : '+str(cubetoread[0]),'yellow')

	weights        = np.arange(0,len(fits),1)
	weights        = np.asarray(Lfir)

	stk_ofn_prfx   = cat_parent + '-' + str(element[2]) + '-' +str(element[1])

	Stack_Res     = Cube_Stack(cubetoread,stk_ofn_prfx,weights,
					sig_clp     = False,sufix=element[0],freq_obs_f=restframe_frequency,
					stack_lite  = stack_light,
					cp_bs_hdrs  = True,
					stt_var     = True,
					spc_wdt_dir = element[0],
					stt_mst_tbl = Cat_Ipt_Tbl,stt_hdr=element[2])
	#############################Add Headers to Stack Results##############################
	name = cat_parent + '-' + str(element[1])
	name = cat_parent + '-' + str(element[2]) + '-' +str(element[1])
	bs_func = ''
	sufix = element[0]
	spc_dir_dst = stk_dir_res

	print
	print "\n".join([file for file in Stack_Res])
	print

	##############################Add Headers to Stack Results##############################
```
This will generate the following fits files in the results directory (```~/Example/Stack_Results-13CO-3D/STACKS/250/```):

```
 - CII_HATLAS-RDS-0-stk-avg-250kms.fits 
 - CII_HATLAS-RDS-0-stk-med-250kms.fits 
 - CII_HATLAS-RDS-0-stk-sum-250kms.fits
```

If ```stack_lite = False``` additional compsoiite spectra will be gnerated:

```
 - CII_HATLAS-RDS-0-stk-3sh-250kms.fits   
 - CII_HATLAS-RDS-0-stk-hsw-250kms.fits   
 - CII_HATLAS-RDS-0-stk-1sh-250kms.fits   
 - CII_HATLAS-RDS-0-stk-3sl-250kms.fits   
 - CII_HATLAS-RDS-0-stk-suw-250kms.fits   
 - CII_HATLAS-RDS-0-stk-1sl-250kms.fits   
 - CII_HATLAS-RDS-0-stk-p25-250kms.fits   
 - CII_HATLAS-RDS-0-stk-wsu-250kms.fits   
 - CII_HATLAS-RDS-0-stk-2sh-250kms.fits   
 - CII_HATLAS-RDS-0-stk-avw-250kms.fits   
 - CII_HATLAS-RDS-0-stk-p75-250kms.fits   
 - CII_HATLAS-RDS-0-stk-2sl-250kms.fits   
 - CII_HATLAS-RDS-0-stk-hst-250kms.fits   
 - CII_HATLAS-RDS-0-stk-std-250kms.fits   
```


###### "Stamps"
Stamps will be created:

```python
		Slices_Files = Cube_Spatial_Extract_Circular(cube2bplot,
													X0_F,Y0_F,
													mask_radi_px_in,mask_radi_as_in,
													mask_radi_px_ot,mask_radi_as_ot,
													mask_radi_px_ms,mask_radi_as_ms,
													x_ref=X0_F_0,y_ref=Y0_F_0,verbose=True,plt_slcs=True,frq_r=restframe_frequency, prefix=prefix_line,
													Splt_Hdr_Cmt_cp=element[2],
													dest_dir_stp = stp_dir_res)
		Plot_Cube_Slices(Slices_Files[0],Slices_Files[1],Slices_Files[2],
						Slices_Files[3],Slices_Files[4],Slices_Files[5],
						frq_r=restframe_frequency, prefix=prefix_line,dest_dir_plt=stm_dir_plt)
```
![Alt text](./Figures/Cube-Slices1.jpg?raw=true "Stacked spectra computed COSMOS field.")
![Alt text](./Figures/Cube-Slices2.jpg?raw=true "Stacked spectra computed COSMOS field.")

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

