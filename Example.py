import sys, os
import numpy as np
from numpy import mean,median
from progressbar import *
from termcolor import colored

from Fnc_Stk_Utl import *
#from Fnc_Stk_Fts import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Tbl import *
from Fnc_Stk_Spc import *
from Fnc_Stk_Stk import *
from Fnc_Stk_Plt import *


os.system('clear')


##################################STACKS#############################################
####################################1################################################
for element in itlpd(channel_width,sbsmn,sbsms):

	nchan         = (2*subcube_width /element[0])+1
	slice_nmb     = int(np.ceil(((nchan-1)/2)))
	#Results
	stk_hme_dir  = home + 'Stack_Results-'+ line +'-3D/'
	img_dir_res  = stk_hme_dir + 'IMAGES/'  #+ str(element[0]) +'/'
	stp_dir_res  = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res  = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_res  = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stm_dir_plt  = plt_dir_res + 'STAMPS/' 
	ana_dir_plt  = plt_dir_res + 'ANALYSIS/'
	mcm_dir_plt  = plt_dir_res + 'MCMC/'
	res_dir_plt  = plt_dir_res + 'RESULTS/' 

	stk_dir_res  = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'
	DIR_CAT_IPT  = [cats_dir]
	DIR_SPC_IPT  = [img_dir]
	DIR_RES      = [stk_hme_dir,stp_dir_res,tbl_dir_res,plt_dir_res,stk_dir_res,stm_dir_plt,ana_dir_plt,mcm_dir_plt,res_dir_plt]
	cat_ipt_tbl  =   cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) 
	Check_directories(cat_ipt_tbl,cat_parent,DIR_RES=DIR_RES)

	cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

	print (colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan'))
	print

	Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
	fits          = Cat_Ipt_Tbl[2]
	delta_nu      = Cat_Ipt_Tbl[4]
	z             = Cat_Ipt_Tbl[8]  
	Lfir          = Cat_Ipt_Tbl[11] 
	nu            = Cat_Ipt_Tbl[13] 
	vel           = Cat_Ipt_Tbl[14] 
	num_obj       = len(Cat_Ipt_Tbl[0])

	cubetoread = [(img_dir_res + image_fits+ '.'+str(element[0]) +'kms.fits') for image_fits in (fits)]
	print (colored('Reading files as : '+str(cubetoread[0]),'yellow'))

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
	print("\n".join([file for file in Stack_Res]))
	print

	##############################Add Headers to Stack Results##############################
#########################################STACKS#############################################
###########################################1################################################


###########################################2################################################
##########################################FIT###############################################
######################################EXTRACT-SQUARE-REGION###################################
for element in itlpd(channel_width,sbsmn,sbsms,func):
	#Results
	stk_hme_dir = home + 'Stack_Results-'+ line +'-3D/'
	#img_dir_res = stk_hme_dir + 'IMAGES/'  + str(element[0])
	stp_dir_res = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_res = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stm_dir_plt = plt_dir_res + 'STAMPS/'  
	ana_dir_plt = plt_dir_res + 'ANALYSIS/'
	mcm_dir_plt = plt_dir_res + 'MCMC/'    
	res_dir_plt = plt_dir_res + 'RESULTS/' 
	stk_dir_res = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'

	stk_hme_dir_ref = home + 'Stack_Results-'+ '12CO' +'-3D/'
	stp_dir_res_ref = stk_hme_dir_ref + 'STAMPS/'  + str(element[0]) +'/'

	DIR_CAT_IPT = [cats_dir]
	DIR_SPC_IPT = [img_dir]
	DIR_RES     = [
					stk_hme_dir,stp_dir_res,tbl_dir_res,
					plt_dir_res,ana_dir_plt,mcm_dir_plt,res_dir_plt,stm_dir_plt,
					stk_dir_res
				]

	#Check_directories(cat_ipt_tbl,cat_parent,DIR_RES=DIR_RES)

	cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt

	print (colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan'))
	print



	nchan             = (2*subcube_width /element[0])+1
	slice_nmb         = int(np.ceil(((nchan-1)/2)))#-1 #19

	cube2bExt         = stk_dir_res+ 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-'+ str(element[3])+'-'+str(element[0])+'kms.fits'
	scale_deg         = Header_Get(cube2bExt,'CDELT2')
	scale_arcsec      = scale_deg*3600#0.00027777778

	X_center,Y_Center = 128,128

	########################################MASK-CIRCULAR-REGION##################################
	for mask_radi_as_ms,mask_radi_as_in,mask_radi_as_ot in zip(apertures_measu,apertures_inner,apertures_outer):
		cube2bplot       = stk_dir_res + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3]+'-'+str(element[0])+'kms.fits'
		scale_deg        = Header_Get(cube2bplot,'CDELT2')
		scale_arcsec     = scale_deg*3600
		mask_radi_as_ms  = mask_radi_as_ms * 1
		mask_radi_px_ms  = mask_radi_as_ms / scale_arcsec
		mask_radi_as_in  = mask_radi_as_in * 1
		mask_radi_px_in  = mask_radi_as_in / scale_arcsec
		mask_radi_as_ot  = mask_radi_as_ot * 1
		mask_radi_px_ot  = mask_radi_as_ot / scale_arcsec

		X0_F,Y0_F        = 128,128
		X0_F_0,Y0_F_0    = 128,128
		print
		print (cube2bplot)
		print ('Masking circular aperture.')
		print ('Center           : ',X0_F_0,Y0_F_0)
		print ('Radii int[arcsec]: ',mask_radi_as_in)
		print ('Radii int[pixels]: ',mask_radi_px_in)
		print ('Radii out[arcsec]: ',mask_radi_as_ot)
		print ('Radii out[pixels]: ',mask_radi_px_ot)
		print ('Radii msr[arcsec]: ',mask_radi_as_ms)
		print ('Radii msr[pixels]: ',mask_radi_px_ms)
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
		########################################MASK-CIRCULAR-REGION###################################
		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + str(element[1]) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[1]) + tbl_ext_ipt
		cat_tbl       = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1]) + tbl_ext_ipt
		print (colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan'))
		print

		Cat_Ipt_Tbl   = Table_Read(cat_tbl,tbl_format_ipt)
		fits          = Cat_Ipt_Tbl[2]
		delta_nu      = Cat_Ipt_Tbl[4]
		z             = Cat_Ipt_Tbl[8]  #[5]
		Lfir          = Cat_Ipt_Tbl[11] #[6]
		nu            = Cat_Ipt_Tbl[13] #[14]
		vel           = Cat_Ipt_Tbl[14] #[15]
		num_obj       = len(Cat_Ipt_Tbl[0])

		z_sample_avg  = np.mean(z)
		z_sample_med  = np.median(z)
		print
		print ('Redshift (avg): ',z_sample_avg)
		print ('Redshift (med): ',z_sample_med)
		print ('subcube_width : ',subcube_width)
		print
		#########################################PLOT-FREQ PROFILE ON REGION############################
		cube2bplot1  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ot) + 'as_msk_ot.fits'
		cube2bplot2  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_in) + 'as_msk_in.fits'
		cube2bplot3  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ms) + 'as_msk_ms.fits'

		cube2bplot4  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ot) + 'as_dta_ot.fits'
		cube2bplot5  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_in) + 'as_dta_in.fits'
		cube2bplot6  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ms) + 'as_dta_ms.fits'

		cube2bplot3_ref = stp_dir_res_ref + '12CO-' + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ms) + 'as_msk_ms.fits'
		cube2bplot6_ref = stp_dir_res_ref + '12CO-' + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) + '-stk-' + element[3] + '-' + str(element[0]) + 'kms-crc-' + str(mask_radi_as_ms) + 'as_dta_ms.fits'

		fit_1D_Gaussian(cube2bplot1,verbose=True,amplitude=0.001,mean=-60,stddev=element[0]/2. * np.sqrt(2. * np.log(2.)),slc_nmb=slice_nmb,max_rng=True,cubewdthv=element[0],rst_frq=restframe_frequency,frq_r=restframe_frequency,prefix=prefix_line,dest_dir_plt = ana_dir_plt)
		fit_1D_Gaussian(cube2bplot2,verbose=True,amplitude=0.001,mean=-60,stddev=element[0]/2. * np.sqrt(2. * np.log(2.)),slc_nmb=slice_nmb,max_rng=True,cubewdthv=element[0],rst_frq=restframe_frequency,frq_r=restframe_frequency,prefix=prefix_line,dest_dir_plt = ana_dir_plt)
		fit_1D_Gaussian(cube2bplot3,verbose=True,amplitude=0.001,mean=-60,stddev=element[0]/2. * np.sqrt(2. * np.log(2.)),slc_nmb=slice_nmb,max_rng=True,cubewdthv=element[0],rst_frq=restframe_frequency,frq_r=restframe_frequency,prefix=prefix_line,dest_dir_plt = ana_dir_plt)
		fit_1D_Gaussian(cube2bplot6,verbose=True,amplitude=0.001,mean=-60,stddev=element[0]/2. * np.sqrt(2. * np.log(2.)),slc_nmb=slice_nmb,max_rng=True,cubewdthv=element[0],rst_frq=restframe_frequency,frq_r=restframe_frequency,prefix=prefix_line,dest_dir_plt = ana_dir_plt)

		Cube_fit_1D_Gaussian(cube2bplot3,
							Cube2bPlot_1D_Err  = cube2bplot1        ,verbose = True   ,
							amplitude          = 0.001              ,mean    = -60    ,stddev    = element[0]/2. * np.sqrt(2. * np.log(2.)),
							slc_nmb            = slice_nmb          ,max_rng = True   ,cubewdthv = element[0]  ,
							rst_frq            = restframe_frequency,frq_r   = restframe_frequency,
							fit_max_1d         = False              ,							
							fit_type           = 'scipy',
							prefix             = line+'-'           ,
							dest_dir_plt       = ana_dir_plt)
		Cube_fit_1D_Gaussian(cube2bplot6,
							Cube2bPlot_1D_Err  = cube2bplot4        ,verbose = True   ,
							amplitude          = 0.001              ,mean    = -60    ,stddev    = element[0]/2. * np.sqrt(2. * np.log(2.)),
							slc_nmb            = slice_nmb          ,max_rng = True   ,cubewdthv = element[0]  ,
							rst_frq            = restframe_frequency,frq_r   = restframe_frequency,
							fit_max_1d         = False              ,							
							fit_type           = 'scipy',
							prefix             = line+'-'           ,
							dest_dir_plt       = ana_dir_plt)
		#########################################PLOT-FREQ PROFILE ON REGION############################
		nchan_ctr = (2*(2000/element[0])+1)
		n_ctr     = (nchan_ctr-1)/2
		print (nchan_ctr)
		print (n_ctr)

		slice_nmb = 8#int(Header_Get(cube2bplot3,'MAX_SNA'))
		slice_nmb = int(Header_Get(cube2bplot3,'MAX_SNA'))
		slice_nmb = n_ctr
		print (slice_nmb)

		#########################################PLOT-STAT CUBE#########################################
		###########################################ON-SLC-REG###########################################
		for function in ['sum','med','avg']: #'sum','med','avg'
			pass
			Plot_Cube_2D(cube2bplot3,
						slc_nmb      = slice_nmb,verbose=True,clp_fnc = str(function),
						redshift     = z_sample_med,rst_frq=restframe_frequency,
						x_ref        = X0_F_0  ,y_ref=Y0_F_0,
						ap_size      = mask_radi_as_ms,
						frq_r        = restframe_frequency,prefix=prefix_line,
						dest_dir_plt = stm_dir_plt,
						dest_dir_clp = stp_dir_res)
		#########################################PLOT-STAT CUBE#########################################
		###########################################ON-SLC-REG###########################################

		#########################################PLOT-STAT CUBE#########################################
		###########################################ON-ORG-CBE###########################################
		Header_Copy(cube2bplot,cube2bplot3,'FTS_FWH')
		Header_Copy(cube2bplot,cube2bplot3,'STT_VEL')
		Header_Copy(cube2bplot,cube2bplot3,'MAX_SNS')
		Header_Copy(cube2bplot,cube2bplot3,'STK_NUM')
		for function in ['sum']: 
			pass
			Plot_Cube_2D(cube2bplot,
						slc_nmb      = slice_nmb,verbose=True,clp_fnc = str(function),
						redshift     = z_sample_med,rst_frq=restframe_frequency,
						x_ref        = X0_F_0,y_ref=Y0_F_0,ap_size=mask_radi_as_in,
						frq_r        = restframe_frequency,prefix=prefix_line,
						dest_dir_plt = stm_dir_plt,
						dest_dir_clp = stp_dir_res)
		#########################################PLOT-STAT CUBE#########################################
		###########################################ON-ORG-CBE###########################################

		#########################################CREATE-SLICES##########################################
		##########################################ON-CLP-CBE############################################
		cubeclp2b_stmp       = stp_dir_res + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3]+'-'+str(element[0])+'kms-2DC-sum.fits'

		Header_Copy(cubeclp2b_stmp,cube2bplot3,'FTS_FWH')
		Header_Copy(cubeclp2b_stmp,cube2bplot3,'STT_VEL')
		Header_Copy(cubeclp2b_stmp,cube2bplot3,'MAX_SNS')
		Header_Copy(cubeclp2b_stmp,cube2bplot3,'STK_NUM')

		Slices_Files = Cube_Spatial_Extract_Circular_2D(cubeclp2b_stmp,
														X0_F,Y0_F,
														mask_radi_px_in,mask_radi_as_in,
														mask_radi_px_ot,mask_radi_as_ot,
														mask_radi_px_ms,mask_radi_as_ms,
														x_ref=X0_F_0,y_ref=Y0_F_0,
														verbose=True,
														frq_r=restframe_frequency, prefix=prefix_line,
														Splt_Hdr_Cmt_cp=element[2],
														dest_dir_stp = stp_dir_res)
		#########################################CREATE-SLICES##########################################
		##########################################ON-CLP-CBE############################################

		###############################FIT 2D GAUSSIAN ON COLLAPSED CUBE################################
		for function in ['sum']:
			pass
			Cube_fit_2D_Gaussian_Noise(cube2bplot3,
										slc_nmb      = None               ,clp_fnc     = function           ,
										sgm_fnc      = element[3]         , 
										SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm)      ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm),
										displ_s_f    = True               ,verbose     = True,circular=True ,
										x_ref        = X0_F_0             ,y_ref       = Y0_F_0             ,
										rst_frq      = restframe_frequency,frq_r       = restframe_frequency,
										sgm_wgth_tms = 'slice_1fw',
										dest_dir_plt = ana_dir_plt,
										dest_dir_clp = stp_dir_res,
										ref_wdt_lne  = spc_wdl_ref       ,ref_wdt_fle = cube2bplot3_ref)
			Cube_fit_2D_Gaussian_Noise(cube2bplot6,
										slc_nmb      = None               ,clp_fnc     = function ,
										SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm)      ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm),
										displ_s_f    = True               ,verbose     = True,circular=True,
										x_ref        = X0_F_0             ,y_ref       = Y0_F_0,
										rst_frq      = restframe_frequency,frq_r       = restframe_frequency,
										sgm_wgth_tms = 'slice_1fw',
										dest_dir_plt = ana_dir_plt,
										dest_dir_clp = stp_dir_res,
										ref_wdt_lne  = spc_wdl_ref       ,ref_wdt_fle = cube2bplot6_ref)
			cube2bplotmskerr = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-crc-'+str(mask_radi_as_ms) + 'as_msk_ms-2DC-'+function+'-nse.fits'
			cube2bplotdtaerr = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-crc-'+str(mask_radi_as_ms) + 'as_dta_ms-2DC-'+function+'-nse.fits'
			slice_nmb00      = 8
			slice_nmb01      = 8
			slice_nmb00      = int(Header_Get(cube2bplot3,'MAX_SNA'))
			slice_nmb01      = int(Header_Get(cube2bplot6,'MAX_SNA'))
			slice_nmb00      = n_ctr
			slice_nmb01      = n_ctr						
			Cube_fit_2D_Gaussian(cube2bplot3,
								Cube2bFit_Err = cube2bplotdtaerr,
								slc_nmb       = None               ,clp_fnc     = function,
								sgm_fnc       = element[3]         , 
								SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm)       ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm)      ,
								displ_s_f     = True               ,verbose     = True,circular=True ,
								x_ref         = X0_F_0             ,y_ref       = Y0_F_0,
								rst_frq       = restframe_frequency,frq_r       = restframe_frequency,
								sgm_wgth_tms  = 'slice_1fw',
								fit_type      = 'scipy'            ,src_sze_fxd = fixed_size,
								dest_dir_plt  = ana_dir_plt,
								ref_wdt_lne   = spc_wdl_ref        ,ref_wdt_fle = cube2bplot3_ref,
								Splt_Hdr_Cmt_cp=element[2]          ,dest_dir_clp  = stp_dir_res)
			Cube_fit_2D_Gaussian(cube2bplot3,
								Cube2bFit_Err = cube2bplotdtaerr,
								slc_nmb       = slice_nmb00        ,clp_fnc     = function,
								sgm_fnc       = element[3]         ,
								SIGMAX_f2DG  = fwhm2sigma(9*tms_sgm)       ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm)      ,
								displ_s_f     = True               ,verbose     = True,circular=True ,
								x_ref         = X0_F_0             ,y_ref       = Y0_F_0,
								rst_frq       = restframe_frequency,frq_r       = restframe_frequency,
								sgm_wgth_tms  = 'slice_1fw',
								fit_type      = 'scipy'            ,src_sze_fxd = fixed_size,
								dest_dir_plt  = ana_dir_plt,
								ref_wdt_lne   = spc_wdl_ref        ,ref_wdt_fle = cube2bplot3_ref,
								Splt_Hdr_Cmt_cp=element[2]          ,dest_dir_clp  = stp_dir_res)
			Cube_fit_2D_Gaussian(cube2bplot6,
								Cube2bFit_Err = cube2bplotmskerr,				
								slc_nmb       = None               ,clp_fnc     = function ,
								sgm_fnc       = element[3]         ,
								SIGMAX_f2DG   = fwhm2sigma(9*tms_sgm)      ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm)      ,
								displ_s_f     = True               ,verbose     = True,circular=True,
								x_ref         = X0_F_0             ,y_ref       = Y0_F_0,
								rst_frq       = restframe_frequency,frq_r       = restframe_frequency,
								sgm_wgth_tms  = 'slice_1fw',
								fit_type      = 'scipy'            ,src_sze_fxd = fixed_size,
								dest_dir_plt  = ana_dir_plt,
								ref_wdt_lne   = spc_wdl_ref        ,ref_wdt_fle = cube2bplot6_ref,
								Splt_Hdr_Cmt_cp=element[2]          ,dest_dir_clp  = stp_dir_res)
			Cube_fit_2D_Gaussian(cube2bplot6,
								Cube2bFit_Err = cube2bplotmskerr,				
								slc_nmb       = slice_nmb01        ,clp_fnc     = function ,
								sgm_fnc       = element[3]         ,
								SIGMAX_f2DG   = fwhm2sigma(9*tms_sgm)      ,SIGMAY_f2DG = fwhm2sigma(9*tms_sgm)      ,
								displ_s_f     = True               ,verbose     = True,circular=True,
								x_ref         = X0_F_0             ,y_ref       = Y0_F_0,
								rst_frq       = restframe_frequency,frq_r       = restframe_frequency,
								sgm_wgth_tms  = 'slice_1fw',#
								fit_type      = 'scipy'            ,src_sze_fxd = fixed_size,
								dest_dir_plt  = ana_dir_plt,
								ref_wdt_lne   = spc_wdl_ref        ,ref_wdt_fle = cube2bplot6_ref,
								Splt_Hdr_Cmt_cp=element[2]          ,dest_dir_clp  = stp_dir_res)
#quit()
		###############################FIT 2D GAUSSIAN ON COLLAPSED CUBE################################
###########################################2################################################
##########################################FIT###############################################

####################################3################################################
##################################Cube Stat###########################################

prefix_line     = line + '-'
func            = ['med','avg']   #sum med avg

if line == '13CO':
	restframe_frequency      =   110.20137E9           
elif line == '12CO':
	restframe_frequency      =   115.271208E9
elif line == '18CO':
	restframe_frequency      =   109.78217340E9

for element in itlpd(channel_width,sbsmn,sbsms,func):	

	cat_ipt_tbl     = cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(element[2]) + '-' +str(element[1])
	cat_tbl         = cat_ipt_tbl + tbl_ext_ipt#

	print
	print (colored('Info from table: ' + str(cat_tbl) + ' ' + str(tbl_format_ipt),'cyan'))
	print


	Cat_Ipt_Tbl = Table_Read(cat_tbl,tbl_format_ipt)
	Cat_Ipt_Stt = Table_Ipt_Cat_Stats(Cat_Ipt_Tbl,element[2])

	z_sample_avg  = Cat_Ipt_Stt[1][0]
	z_sample_med  = Cat_Ipt_Stt[1][1]
	print
	print ('Redshift (avg): ',z_sample_avg)
	print ('Redshift (med): ',z_sample_med)
	print ('subcube_width : ',subcube_width)

	#Results
	stk_hme_dir = home + 'Stack_Results-'+ line +'-3D/'
	#img_dir_res = stk_hme_dir + 'IMAGES/'  + str(channel_width)
	stp_dir_res = stk_hme_dir + 'STAMPS/'  + str(element[0]) +'/'
	tbl_dir_res = stk_hme_dir + 'TABLES/'  + str(element[0]) +'/'
	plt_dir_res = stk_hme_dir + 'PLOTS/'   + str(element[0]) +'/'
	stk_dir_res = stk_hme_dir + 'STACKS/'  + str(element[0]) +'/'

	DIR_CAT_IPT = [cats_dir]
	DIR_SPC_IPT = [img_dir]
	DIR_RES     = [
				stk_hme_dir,stp_dir_res,tbl_dir_res,
				plt_dir_res,ana_dir_plt,mcm_dir_plt,res_dir_plt,
				stk_dir_res
			]

	Check_directories(cat_ipt_tbl,cat_parent)


	apertures_measu = [15]
	apertures_inner = [10]
	apertures_outer = [20]

	nchan_ctr = (2*(2000/element[0])+1)
	n_ctr     = (nchan_ctr-1)/2

	for mask_radi_as_ms,mask_radi_as_in,mask_radi_as_ot in zip(apertures_measu,apertures_inner,apertures_outer):
		cube2bplot_ot  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-crc-'+str(mask_radi_as_ot) + 'as_msk_ot.fits' 
		cube2bplot_in  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-crc-'+str(mask_radi_as_in) + 'as_msk_in.fits' 
		cube2bplot_ms  = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-crc-'+str(mask_radi_as_ms) + 'as_msk_ms.fits' 
		#
		cube2bplot_clp_ot = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-2DC-sum-crc-'+str(mask_radi_as_ot) + 'as_msk_ot.fits'
		cube2bplot_clp_in = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-2DC-sum-crc-'+str(mask_radi_as_in) + 'as_msk_in.fits'
		cube2bplot_clp_ms = stp_dir_res + prefix_line + 'CII_HATLAS-'+ str(element[2]) + '-' +str(element[1]) +'-stk-'+ element[3] +'-'+str(element[0])+'kms-2DC-sum-crc-'+str(mask_radi_as_ms) + 'as_msk_ms.fits'
		#
		slc_nmb1 = n_ctr
		slc_nmb1 = n_ctr
		Cube_Stat(cube2bplot_in,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency,dest_dir_tbl=tbl_dir_res)
		Cube_Stat(cube2bplot_ot,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency,dest_dir_tbl=tbl_dir_res)
		Cube_Stat(cube2bplot_ms,redshift=z_sample_med,rst_frq=restframe_frequency,slc_nmb=slc_nmb1,cubewdthv=int(element[0]),frq_r=restframe_frequency,dest_dir_tbl=tbl_dir_res)
##########################################Cube Stat##########################################
############################################3################################################

############################################4#################################################
#############################################MCMC#############################################
method          = 3          #1-sum 2-1DGF 3-2DGF-sigma 4-2DGF-fwhm 5-2DGF
error           = 1          #N sigma confidence intervals
plt_scl         = None       #Z, Y, both None
log_lm          = 'error'    #both,value,error, None
func1           = 'avg'
func2           = 'med'
cw              = 250
sbsmn           = [0]
sbsms           = 'RDS_B'
mask_radi_as_ms = 15

ERR_MC_CLC      = False
ERR_MC_ERR_PLT  = True

lit_res         = False   #Plot literature Results

MCMC_generator(iterations_mc,line1,line2,method,error,'RDS_B',sbsmn,spc_wdt_dir=cw,mask_radi_as_ms=mask_radi_as_ms)
quit()
