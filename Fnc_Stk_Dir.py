import sys, os
import numpy as np
from numpy import mean,median
from progressbar import *
from termcolor import colored
from astropy import constants as apct
from astropy import units as u
from itertools import product as itlpd

from astropy.io.fits import getdata as apgtdt
from astropy.wcs import WCS as apwcs
from astropy import cosmology as apcosmo

from spectral_cube import SpectralCube as scspc

from Fnc_Stk_Utl import *

import platform
py_ver = (platform.python_version_tuple())
py_ver_ma = py_ver[0]
py_ver_mc = py_ver[1]
py_ver_mx = py_ver[2]
print
print (colored('Python version: ' + str(py_ver_ma) + '.' + str(py_ver_mc) +'.' +  str(py_ver_mx),'yellow'))
print


home = home = os.path.expanduser("~")+ '/Desktop/Example-VSAT-3D/'
line = '13CO' 

def fwhm2sigma(fwhm):
	return fwhm/(2*np.sqrt(2*np.log(2)))
def sigma2fwhm(sigma):
	return sigma*(2*np.sqrt(2*np.log(2)))

#Catalogue
cat_parent               =   'CII_HATLAS'
CAT_PARENT               =   cat_parent.upper()

#Image Directories
cats_dir                 =   home     + '/Catalogues/'
cat_dir                  =   cats_dir + CAT_PARENT + '/'
dts_dir                  =   home     + line +'DataSets/'             	 #DataSet_directory
ext_dir                  =   home     + line +'Extracted-Im/'         	 #Extracted_directory
fts_dir                  =   ext_dir  + 'FITS/'                       	 #Fits_directory
img_dir                  =   ext_dir  + 'IMAGES/'                     	 #Images_directory

#Frequecies
subcube_width            =   2000                              		  	 #kms-1
channel_width            =   [250]                             		  	 #kms-1

#Input Table
sbsmn                    =   [0]
sbsms                    =   ['RDS_B']						 			#STM,SFR,LCO,sSF,MH2,SFE,SDG,SDS,TDT,RDS
cat_ipt_tbl              =   cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(sbsms) + '-' +str(sbsmn) 
unq_tbl_ipt              =   'yes'
head_un_ipt              =   'ident'									#'spec1d'

#Stacking
stack_light     		 = 	 True 							 			#True: SUM MED AVG False: Histograms and pct (1,2,3-sigma)stacks

#Sigma-Clip
sigma_clipping           =   True                            			# Sigma clipping
sigma_cut                =   0                               			# sigma cut
sigma_cen_fct            =   mean                            			# median, mean
sigma_msk_fill_val       =   np.nan                          			# np.nan, value

#Fitting
func            		=	 ['avg','med']  				 			#sum med avg
apertures_measu 		=	 [15] 							 			#radii in as measure circle
apertures_inner 		=	 [10] 							 			#radii in as inner circle
apertures_outer 		=	 [20] 							 			#radii in as outer circle
prefix_line     		=	 line + '-'
spc_wdl_ref     		=	 False                			 			#Use 12CO/other fit width line for Flux computations cube2bplot6_ref
fixed_size      		=	 True                 			 			#Fix Source Size
tms_sgm         		=	 1                    			 			#SOURCE SIZE BASED ON THE SYN BEAM IN SIGMA TERMS

#MCMC
iterations_mc   		 =   1000   									#50000 #mcmc_itr_nmb
plot_dist_hist  		 =   True  										#Plot MCMC histograms
line1  			         =   '12CO'
line2  			         =   '13CO'

#tables
tbl_format_ipt           =   'csv'                           			#ascii,csv,fits,ascii.fixed_width_two_line
tbl_format_opt           =   'csv'                           			#ascii,csv,fits,ascii.fixed_width_two_line

#Results
stk_hme_dir              =   home + 'Stack_Results-'+ line +'-3D/'
img_dir_res              =   stk_hme_dir + 'IMAGES/' 
stp_dir_res              =   stk_hme_dir + 'STAMPS/' 
tbl_dir_res              =   stk_hme_dir + 'TABLES/' 
plt_dir_tbl              =   tbl_dir_res + 'PLOTS/' 
plt_dir_res              =   stk_hme_dir + 'PLOTS/'
stm_dir_plt              =   plt_dir_res + 'STAMPS/' 
ana_dir_plt              =   plt_dir_res + 'ANALYSIS/'
mcm_dir_plt              =   plt_dir_res + 'MCMC/'
res_dir_plt              =   plt_dir_res + 'RESULTS/' 
stk_dir_res              =   stk_hme_dir + 'STACKS/'

#Output  Tables
ind_tbl                  =   'yes'
grl_tbl                  =   'yes'

grl_tbl_nmB              =   'Prs_Bkg_'+ CAT_PARENT 
grl_tbl_nmF              =   'Prs_Frg_'+ CAT_PARENT 

unq_tbl_opt              =   'yes'
hed_un_opt_F             =   'id_F'                                #header of column for uniqueness
hed_un_opt_B             =   'id_B'                                #header of column for uniqueness
grl_tbl_nmB_U            =   'Prs_Bkg_'+ CAT_PARENT +'_U'
grl_tbl_nmF_U            =   'Prs_Frg_'+ CAT_PARENT +'_U'

verbose                  =   'yes'

if line == '13CO':
	restframe_frequency      =   110.20137E9           
elif line == '12CO':
	restframe_frequency      =   115.271208E9
elif line == '18CO':
	restframe_frequency      =   109.78217340E9
#############################################################################################################################
DIR_CAT_IPT = [cats_dir]
DIR_SPC_IPT = [img_dir]
DIR_RES     = [
				stk_hme_dir,img_dir_res,stp_dir_res,tbl_dir_res,
				plt_dir_res,stk_dir_res,plt_dir_tbl,
				ana_dir_plt,mcm_dir_plt,res_dir_plt]

if tbl_format_ipt == 'ascii' or tbl_format_ipt == 'ascii.fixed_width_two_line':
	tbl_ext_ipt = '.dat'
elif tbl_format_ipt == 'csv':	
	tbl_ext_ipt = '.csv'
elif tbl_format_ipt == 'fits':	
	tbl_ext_ipt = '.fits'

if tbl_format_opt == 'ascii' or tbl_format_opt == 'ascii.fixed_width_two_line':
	tbl_ext_opt = '.dat'
elif tbl_format_opt == 'csv':	
	tbl_ext_opt = '.csv'
elif tbl_format_opt == 'fits':	
	tbl_ext_opt = '.fits'

cat_tbl    = cat_ipt_tbl + tbl_ext_ipt
cat_tbl_U  = cat_ipt_tbl + tbl_ext_opt

sofl = apct.c
sofl = sofl.to(u.km/u.s)

cosmo_H0 = 70
cosmo_omegaM = 0.3
cosmo = apcosmo.FlatLambdaCDM(cosmo_H0,cosmo_omegaM)
#############################################################################################################################	
def Check_directories(cat_tbl_chk,cat_chk,*args, **kwargs):
	DIR_RES = kwargs.get('DIR_RES',[
				stk_hme_dir,img_dir_res,stp_dir_res,tbl_dir_res,
				plt_dir_res,stk_dir_res,plt_dir_tbl,
				ana_dir_plt,mcm_dir_plt,res_dir_plt])
	print
	print ('Checking input directories of the catalogue : ',cat_chk)
	if os.path.exists(str(cat_tbl_chk)+str(tbl_ext_ipt))==True :
		print
		print ('Catalogue table exists             : ', str(cat_tbl_chk)+str(tbl_ext_ipt))
		print ('Spectra directories exists         : ', str(DIR_SPC_IPT[0]))
		print
		print ('Checking Result Directories.')
		print

		for tree in DIR_RES:
			if os.path.isdir(tree)==True:
				pass
				print ('Directory exists: ', tree)
			elif os.path.isdir(tree)==False:
				print ('Directory does not exist, creating it: ', tree)
				os.makedirs(tree)
	elif os.path.exists(str(cat_tbl_chk)+str(tbl_ext_ipt))==False :
		print
		print ('Some of the directories does not exist.')
		print ('Check input directories. ')
		print (str(cat_tbl_chk)+str(tbl_ext_ipt))
		print( DIR_SPC_IPT[0])
#############################################################################################################################
