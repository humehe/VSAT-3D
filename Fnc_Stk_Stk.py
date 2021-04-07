import bottleneck as bn
from astropy import stats as apsts

import scipy.integrate as integrate

from Fnc_Stk_Dir import *
from Fnc_Stk_Spc import *
from Fnc_Stk_Tbl import *

####Fnc_Stk_Stk###
def Cube_Stack(Cubes2bStacked,name,wght_img_2bstack,sig_clp,*args, **kwargs):
	wrt_fits         = kwargs.get('wrt_fits'       ,True)
	pst_msk          = kwargs.get('pst_msk'        ,False)
	pst_smt          = kwargs.get('pst_smt'        ,False)
	pst_cnt          = kwargs.get('pst_cnt'        ,False)
	stack_ext        = kwargs.get('stack_ext'      ,None)
	new_CRVAL1_head  = kwargs.get('new_CRVAL1_head',None)
	new_CDELT1_head  = kwargs.get('new_CDELT1_head',None)
	smt_spc_pst      = kwargs.get('smt_spc_pst'    ,False)
	smooth_shape     = kwargs.get('smooth_shape'   ,'gaussian')
	wght_type        = kwargs.get('wght_type'      ,None)
	
	sufix            = kwargs.get('sufix'          ,'')
	freq_obs_f       = kwargs.get('freq_obs_f'     ,99999)

	stack_lite       = kwargs.get('stack_lite'     ,True)

	spc_wdt_dir      = kwargs.get('spc_wdt_dir'    ,500)

	cp_bs_hdrs       = kwargs.get('cp_bs_hdrs'     ,False)

	stt_var			 = kwargs.get('stt_var',False)
	stt_mst_tbl		 = kwargs.get('stt_mst_tbl',None)
	stt_hdr			 = kwargs.get('stt_hdr',None)

	[Cube_Mask_Nan(img)     for img in Cubes2bStacked]

	img_2bstack    = [apgtdt(img,memmap=False) for img in Cubes2bStacked]
	wcs            = kwargs.get('wcs'            ,apwcs(Cubes2bStacked[0]))
	try:
		wcs       = wcs.dropaxis(3)
	except IndexError:
		pass

	print
	print ('Number of galaxies to be stacked (histogram): ',len(img_2bstack))

	if sig_clp == True:
		img_flt       = apsts.sigma_clip(img_2bstack,sigma=sigma_cut,axis=0,iters=None,cenfunc=sigma_cen_fct, copy=True)

		print
		print (colored('Sigma-clipping for stacking!','yellow'))
		print (colored('Sigma Cut                    : ' + str(sigma_cut),'yellow'))
		print (colored('Central function             : ' + str(sigma_cen_fct), 'yellow'))
		print (colored('Central Value for clipping   : ' + str(sigma_cen_fct),'yellow'))

		img_flt.set_fill_value(sigma_msk_fill_val)
		img_flt_filled = img_flt.filled()
		img_stat       = img_flt_filled
	elif sig_clp == False:
		img_stat   = img_2bstack

	wght_img_copy = wght_img_2bstack
	wght_img_stat = wght_img_2bstack
	wght_img_stat = np.asarray(wght_img_stat)
	img_staw      = []

	img_stat_smw_f = []

	[img_staw.append(np.asarray(img_stat)[j]*np.asarray(wght_img_stat)[j]) for j in range(len(wght_img_stat))]
	img_staw      = np.asarray(img_staw)
	[img_stat_smw_f.append(np.divide(np.asarray(img_staw)[j],np.asarray(img_stat)[j])) for j in range(len(img_stat))]

	print
	print (colored('Original shape                                               : '+str(np.asarray(img_stat).shape),'cyan'))
	img_stat = np.squeeze(img_stat) 
	N,F,l,m = np.asarray(img_stat).shape 
	print (colored('Squeezed useless extra dimensions                            : '+str(np.asarray(img_stat).shape),'cyan'))
	print (colored('Dimension Numer of Cubes, Number of channels, X size, Y size : '+str(N)+', '+str(F)+', '+str(l)+', '+str(m),'cyan'))
	print

	img_res_sum = scspc(data=bn.nansum(np.array(img_stat)             , axis=0)  ,wcs=wcs) 
	img_res_avg = scspc(data=bn.nanmean(np.array(img_stat)            , axis=0)  ,wcs=wcs) 
	img_res_med = scspc(data=bn.nanmedian(np.array(img_stat)          , axis=0)  ,wcs=wcs) 

	print
	print (colored('Sum, Mean, Median : Stacked data cubes OK','yellow'))
	print 

	if stack_lite == False:

		img_stat_hst_f = []
		img_stat_hsw_f = []
		#BEGINS HISTO
		widgets = ['Computing histogram for Stacks: ', Percentage(), ' ', Bar(marker='*',left='[',right=']'),
				   ' ', ETA(), ' ', FileTransferSpeed()] 

		pbar = ProgressBar(widgets=widgets, maxval=F)
		pbar.start()

		for freq in range(F):
			pbar.update(freq)
			FREQ = np.asarray(img_stat)[:,freq,:,:]
			img_stat_hst_y = []
			img_stat_hsw_y = []
			for y_dim in range(l):
				Y_ROW = np.asarray(img_stat)[:,freq,y_dim,:]
				Transpose  = np.asarray(Y_ROW).T
				Transposw  = np.asarray(Y_ROW).T
				img_stat_hst_x = []
				img_stat_hsw_x = []
				for x_dim in range(len(Transpose)):
					if np.isnan(sigma_msk_fill_val) == True:
						non_msk_num = int(np.count_nonzero(~np.isnan(Transpose[x_dim])))
						msk_num     = int(np.count_nonzero(np.isnan(Transpose[x_dim])))
						img_stat_hst_x.append(float(non_msk_num))

						non_msk_num_wghts = int(np.count_nonzero(~np.isnan(Transposw[x_dim])))
						msk_num_wghts     = int(np.count_nonzero(np.isnan(Transposw[x_dim])))
						img_stat_hsw_x.append(float(non_msk_num_wghts))

					elif np.isnan(sigma_msk_fill_val) == False:
						pass
						non_msk_num = int(np.count_nonzero(Transpose[x_dim]!=sigma_msk_fill_val))
						img_stat_hst_x.append(float(non_msk_num))

						non_msk_num_wghts = int(np.count_nonzero(Transposw[x_dim]!=sigma_msk_fill_val))
						img_stat_hsw_x.append(float(non_msk_num_wghts))
					else:
						pass
				
				img_stat_hst_x = np.reshape(img_stat_hst_x,(m))
				img_stat_hsw_x = np.reshape(img_stat_hsw_x,(m))

				img_stat_hst_y.append(img_stat_hst_x)
				img_stat_hsw_y.append(img_stat_hsw_x)

			img_stat_hst_f.append(img_stat_hst_y)
			img_stat_hsw_f.append(img_stat_hsw_y)
		pbar.finish()
		#ENDS HISTO

		img_sts_hst = scspc(data=np.asarray(img_stat_hst_f)                          ,wcs=wcs) 
		img_res_std = scspc(data=bn.nanstd(np.array(img_stat)             , axis=0)  ,wcs=wcs)

		print
		print (colored('Histogram, Std: Stacked data cubes OK','yellow'))
		print 

		img_res_suw_pre = np.asarray(bn.nansum(np.array(img_staw)                , axis=0))
		img_sts_wsu_pre = np.asarray(bn.nansum(np.array(img_stat_smw_f)          , axis=0))

		img_sts_wsu_pre = np.squeeze(img_sts_wsu_pre)
		img_res_suw_pre = np.squeeze(img_res_suw_pre)

		print
		print (colored('Weights Sum Weighted Sum pre computations: OK','yellow'))
		print 

		img_sts_hsw = scspc(data=np.asarray(img_stat_hsw_f)                          ,wcs=wcs)
		img_sts_wsu = scspc(data=img_sts_wsu_pre                                     ,wcs=wcs)
		img_res_suw = scspc(data=img_res_suw_pre                                     ,wcs=wcs)
		img_res_avw = scspc(data=img_res_suw_pre.astype(float)/img_sts_wsu_pre.astype(float) ,wcs=wcs)
		
		print
		print (colored('SW Histogram, Sum of weights, Weighted Sum: Stacked data cubes OK','yellow'))
		print

		img_res_1sl = scspc(data=np.nanpercentile(np.array(img_stat), 15.9, axis=0)  ,wcs=wcs)	
		img_res_1sh = scspc(data=np.nanpercentile(np.array(img_stat), 84.1, axis=0)  ,wcs=wcs)
		img_res_2sl = scspc(data=np.nanpercentile(np.array(img_stat), 2.30, axis=0)  ,wcs=wcs)
		img_res_2sh = scspc(data=np.nanpercentile(np.array(img_stat), 97.7, axis=0)  ,wcs=wcs)
		img_res_3sl = scspc(data=np.nanpercentile(np.array(img_stat), 0.20, axis=0)  ,wcs=wcs)
		img_res_3sh = scspc(data=np.nanpercentile(np.array(img_stat), 99.8, axis=0)  ,wcs=wcs)
		img_res_p25 = scspc(data=np.nanpercentile(np.array(img_stat), 25.0, axis=0)  ,wcs=wcs)
		img_res_p75 = scspc(data=np.nanpercentile(np.array(img_stat), 75.0, axis=0)  ,wcs=wcs)

		print ('Stacked images through : sum, mean, median, and percentiles: ')
		print ('17., 83.0, (1 sigma)')
		print ('2.5, 97.5, (2 sigma)')
		print ('0.5, 99.5, (3 sigma)')
		print ('25., 75.0, (interquantile)')
		print		
		print (colored('Percentiles: Stacked data cubes OK','yellow'))
		print 
	elif stack_lite == True:
		pass

	bs_func = kwargs.get('bs_func','')

	if wrt_fits==True:
		if  '-BS-' in name:
			print (colored(name,'yellow'))
			spc_dir_dst = str_bst_stk + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass
		elif  '-BS_MST' in name:
			print (colored(name,'yellow'))
			spc_dir_dst = stt_bst_stk + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass
		else:
			spc_dir_dst = stk_dir_res + str(spc_wdt_dir) +'/'
			if os.path.exists(spc_dir_dst)==False:
				print
				print (colored('Stacked width directory does not exist!','yellow'))
				print (colored('Creating it!','yellow'))
				print
				os.makedirs(spc_dir_dst)
			else:
				pass

		spec_file_sum_ofn = spc_dir_dst + str(name) + bs_func + '-stk-sum-' + str(sufix) + 'kms.fits'
		spec_file_avg_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avg-' + str(sufix) + 'kms.fits'
		spec_file_med_ofn = spc_dir_dst + str(name) + bs_func + '-stk-med-' + str(sufix) + 'kms.fits'
		spec_file_hst_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hst-' + str(sufix) + 'kms.fits'
		spec_file_std_ofn = spc_dir_dst + str(name) + bs_func + '-stk-std-' + str(sufix) + 'kms.fits'
		spec_file_p25_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p25-' + str(sufix) + 'kms.fits'
		spec_file_p75_ofn = spc_dir_dst + str(name) + bs_func + '-stk-p75-' + str(sufix) + 'kms.fits'
		spec_file_1sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sl-' + str(sufix) + 'kms.fits'
		spec_file_1sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-1sh-' + str(sufix) + 'kms.fits'
		spec_file_2sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sl-' + str(sufix) + 'kms.fits'
		spec_file_2sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-2sh-' + str(sufix) + 'kms.fits'
		spec_file_3sl_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sl-' + str(sufix) + 'kms.fits'
		spec_file_3sh_ofn = spc_dir_dst + str(name) + bs_func + '-stk-3sh-' + str(sufix) + 'kms.fits'
		
		spec_file_hsw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-hsw-' + str(sufix) + 'kms.fits'
		spec_file_wsu_ofn = spc_dir_dst + str(name) + bs_func + '-stk-wsu-' + str(sufix) + 'kms.fits'
		spec_file_suw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-suw-' + str(sufix) + 'kms.fits'
		spec_file_avw_ofn = spc_dir_dst + str(name) + bs_func + '-stk-avw-' + str(sufix) + 'kms.fits'

		spec_file_sum     = img_res_sum.write(spec_file_sum_ofn,overwrite=True)
		spec_file_avg     = img_res_avg.write(spec_file_avg_ofn,overwrite=True)
		spec_file_med     = img_res_med.write(spec_file_med_ofn,overwrite=True)

		Cube_Freq2VelAxis(spec_file_sum_ofn)
		Cube_Freq2VelAxis(spec_file_avg_ofn)
		Cube_Freq2VelAxis(spec_file_med_ofn)

		if stack_lite == False:
			spec_file_hst     = img_sts_hst.write(spec_file_hst_ofn,overwrite=True)
			spec_file_std     = img_res_std.write(spec_file_std_ofn,overwrite=True)

			spec_file_p25     = img_res_p25.write(spec_file_p25_ofn,overwrite=True)
			spec_file_p75     = img_res_p75.write(spec_file_p75_ofn,overwrite=True)
			spec_file_1sl     = img_res_1sl.write(spec_file_1sl_ofn,overwrite=True)
			spec_file_1sh     = img_res_1sh.write(spec_file_1sh_ofn,overwrite=True)
			spec_file_2sl     = img_res_2sl.write(spec_file_2sl_ofn,overwrite=True)
			spec_file_2sh     = img_res_2sh.write(spec_file_2sh_ofn,overwrite=True)
			spec_file_3sl     = img_res_3sl.write(spec_file_3sl_ofn,overwrite=True)
			spec_file_3sh     = img_res_3sh.write(spec_file_3sh_ofn,overwrite=True)

			spec_file_hsw     = img_sts_hsw.write(spec_file_hsw_ofn,overwrite=True)
			spec_file_wsu     = img_sts_wsu.write(spec_file_wsu_ofn,overwrite=True)
			spec_file_suw     = img_res_suw.write(spec_file_suw_ofn,overwrite=True)

			spec_file_avw     = img_res_avw.write(spec_file_avw_ofn,overwrite=True)


			Cube_Freq2VelAxis(spec_file_hst_ofn)
			Cube_Freq2VelAxis(spec_file_std_ofn)

			Cube_Freq2VelAxis(spec_file_p25_ofn)
			Cube_Freq2VelAxis(spec_file_p75_ofn)
			Cube_Freq2VelAxis(spec_file_1sl_ofn)
			Cube_Freq2VelAxis(spec_file_1sh_ofn)
			Cube_Freq2VelAxis(spec_file_2sl_ofn)
			Cube_Freq2VelAxis(spec_file_2sh_ofn)
			Cube_Freq2VelAxis(spec_file_3sl_ofn)
			Cube_Freq2VelAxis(spec_file_3sh_ofn)
			Cube_Freq2VelAxis(spec_file_hsw_ofn)
			Cube_Freq2VelAxis(spec_file_wsu_ofn)
			Cube_Freq2VelAxis(spec_file_suw_ofn)
			Cube_Freq2VelAxis(spec_file_avw_ofn)

			OPT_STCK_FLS = [spec_file_sum_ofn,spec_file_avg_ofn,spec_file_med_ofn,spec_file_hst_ofn,
			spec_file_std_ofn,
			spec_file_p25_ofn,spec_file_p75_ofn,
			spec_file_1sl_ofn,spec_file_1sh_ofn,
			spec_file_2sl_ofn,spec_file_2sh_ofn,
			spec_file_3sl_ofn,spec_file_3sh_ofn,
			spec_file_hsw_ofn,spec_file_wsu_ofn,spec_file_suw_ofn,spec_file_avw_ofn]
		elif stack_lite == True:
			OPT_STCK_FLS = [spec_file_sum_ofn,spec_file_avg_ofn,spec_file_med_ofn]
		[Header_Updt(spec_sts_res,'STK_NUM' ,len(img_2bstack), header_comment = 'Number of galaxies used for Stack') for spec_sts_res in OPT_STCK_FLS]
	else:
		pass


	print ('Imaged Stacked files names: ')
	print
	print (colored(spec_file_sum_ofn,'cyan'))
	print (colored(spec_file_avg_ofn,'cyan'))
	print (colored(spec_file_med_ofn,'cyan'))
	if stack_lite == False:
		print (colored(spec_file_hst_ofn,'cyan'))
		print (colored(spec_file_std_ofn,'cyan'))
		print (colored(spec_file_p25_ofn,'cyan'))
		print (colored(spec_file_p75_ofn,'cyan'))
		print (colored(spec_file_1sl_ofn,'cyan'))
		print (colored(spec_file_1sh_ofn,'cyan'))
		print (colored(spec_file_2sl_ofn,'cyan'))
		print (colored(spec_file_2sh_ofn,'cyan'))
		print (colored(spec_file_3sl_ofn,'cyan'))
		print (colored(spec_file_3sh_ofn,'cyan'))
		

		print (colored(spec_file_hsw_ofn,'yellow'))
		print (colored(spec_file_avw_ofn,'yellow'))
		print (colored(spec_file_suw_ofn,'yellow'))
	elif stack_lite == True:
		pass

	if stack_lite == True:
		FNL_SPEC_RES = [spec_file_med,spec_file_avg,spec_file_sum]
	elif stack_lite == False:
		FNL_SPEC_RES = [
					spec_file_med,spec_file_avg,spec_file_sum,spec_file_std,
					spec_file_hst,
					spec_file_1sl,spec_file_1sh,
					spec_file_2sl,spec_file_2sh,
					spec_file_3sl,spec_file_3sh,
					spec_file_p25,spec_file_p75,
					spec_file_hsw,spec_file_wsu,spec_file_suw,spec_file_avw]
	if cp_bs_hdrs == True:
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BSCALE')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BZERO')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BMAJ')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BMIN')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BPA')      for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BTYPE')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'EQUINOX')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'RADESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'BUNIT')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'RADESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'LONPOLE')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'LATPOLE')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC1_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC2_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC3_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC1_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC2_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC3_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC1_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC2_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PC3_3')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CTYPE1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRVAL1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CDELT1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRPIX1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CUNIT1')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CTYPE2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRVAL2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CDELT2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRPIX2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CUNIT2')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CTYPE3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRVAL3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CDELT3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CRPIX3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'CUNIT3')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PV2_1')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'PV2_2')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'RESTFRQ')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'SPECSYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'ALTRVAL')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'ALTRPIX')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'VELREF')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'TELESCOP') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSERVER') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'DATE-OBS') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'TIMESYS')  for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSRA')    for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSDEC')   for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSGEO-X') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSGEO-Y') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'OBSGEO-Z') for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'DATE')     for stk_res_flr in OPT_STCK_FLS]
		[Header_Copy(stk_res_flr,Cubes2bStacked[0],'ORIGIN')   for stk_res_flr in OPT_STCK_FLS]
	else:
		pass
	if stt_var == True:
		print
		print (colored('Adding stat to fits headers!','yellow'))
		print
		tbl_sts = Table_Ipt_Cat_Stats(stt_mst_tbl,stt_hdr)
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][0] ,tbl_sts[1][0] ,header_comment='Redshift Average')                         for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][1] ,tbl_sts[1][1] ,header_comment='Redshift Median')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][2] ,tbl_sts[1][2] ,header_comment='Redshift 1 sgm lw lmt 15.9 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][3] ,tbl_sts[1][3] ,header_comment='Redshift 1 sgm hg lmt 84.1 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][4] ,tbl_sts[1][4] ,header_comment='Redshift 2 sgm lw lmt 2.30 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][5] ,tbl_sts[1][5] ,header_comment='Redshift 2 sgm hg lmt 97.7 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][6] ,tbl_sts[1][6] ,header_comment='Redshift 3 sgm lw lmt 0.20 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][7] ,tbl_sts[1][7] ,header_comment='Redshift 3 sgm hg lmt 99.8 pct')           for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][8] ,tbl_sts[1][8] ,header_comment='Redshift 25 pct')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][9] ,tbl_sts[1][9] ,header_comment='Redshift 75 pct')                          for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][10],tbl_sts[1][10],header_comment=str(tbl_sts[2]) + ' Average')               for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][11],tbl_sts[1][11],header_comment=str(tbl_sts[2]) + ' Median')                for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][12],tbl_sts[1][12],header_comment=str(tbl_sts[2]) + ' 1 sgm lw lmt 15.9 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][13],tbl_sts[1][13],header_comment=str(tbl_sts[2]) + ' 1 sgm hg lmt 84.1 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][14],tbl_sts[1][14],header_comment=str(tbl_sts[2]) + ' 2 sgm lw lmt 2.30 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][15],tbl_sts[1][15],header_comment=str(tbl_sts[2]) + ' 2 sgm hg lmt 97.7 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][16],tbl_sts[1][16],header_comment=str(tbl_sts[2]) + ' 3 sgm lw lmt 0.20 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][17],tbl_sts[1][17],header_comment=str(tbl_sts[2]) + ' 3 sgm hg lmt 99.8 pct') for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][18],tbl_sts[1][18],header_comment=str(tbl_sts[2]) + ' 25 pct')                for Stacked_Cube in OPT_STCK_FLS]
		[Header_Get_Add(Stacked_Cube,tbl_sts[0][19],tbl_sts[1][19],header_comment=str(tbl_sts[2]) + ' 75 pct')                for Stacked_Cube in OPT_STCK_FLS]		
	else:
		pass

	return OPT_STCK_FLS
####Fnc_Stk_Stk###	