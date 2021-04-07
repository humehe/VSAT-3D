from astropy import table as aptbl
from scipy.constants import physical_constants

from Fnc_Stk_Dir import *
from Fnc_Stk_Fts import *
from Fnc_Stk_Mth import *

####Fnc_Stk_Stt###
def Cube_Stat(Cube2Stt,*args,**kwargs):
	dest_dir  = kwargs.get('dest_dir', None)
	autoaxis  = kwargs.get('autoaxis', False)
	verbose   = kwargs.get('verbose' , False)
	epssave   = kwargs.get('epssave' , False)
	showplot  = kwargs.get('showplot', False)
	clp_fnc   = kwargs.get('clp_fnc' ,'sum')

	z_avg     = kwargs.get('z_avg',Header_Get(Cube2Stt,'STZ_AVG'))
	z_med     = kwargs.get('z_med',Header_Get(Cube2Stt,'STZ_MED'))
	frq_r     = kwargs.get('frq_r',Header_Get(Cube2Stt,'RESTFRQ'))
	z_f2l     = z_med
	cubewdthv = kwargs.get('cubewdthv'      ,1)


	x_ref    = kwargs.get('x_ref',0)
	y_ref    = kwargs.get('y_ref',0)
	ap_size  = kwargs.get('ap_size',0)

	cube_data     = np.asarray(apgtdt(Cube2Stt,memmap=False) )
	freq_num,y_num,x_num = cube_data.shape

	slc_nmb  = kwargs.get('slc_nmb', freq_num/2)


	TABLESTATNAME_1 = (str(Cube2Stt).split('.fits')[0]).split('/')[-1] + '-stt.dat'
	TABLESTATNAME_2 = (str(Cube2Stt).split('.fits')[0]).split('/')[-1] + '-stt.'+ tbl_format_opt
	if dest_dir != None:
		TABLESTATNAME_1 = str(dest_dir)  + '/' + TABLESTATNAME_1
		TABLESTATNAME_2 = str(dest_dir)  + '/' + TABLESTATNAME_2
	elif dest_dir == None:
		TABLESTATNAME_1 = tbl_dir_res    + '/' + TABLESTATNAME_1
		TABLESTATNAME_2 = tbl_dir_res    + '/' + TABLESTATNAME_2
	Message1 = 'Cube Stat on Table: ' + TABLESTATNAME_1 + ' slice number : ' + str(slc_nmb)
	Message2 = 'Cube Stat on Table: ' + TABLESTATNAME_2 + ' slice number : ' + str(slc_nmb)

	cube_data_sgl_slc = cube_data[slc_nmb]
	cube_data_clp_sum = np.asarray(np.nansum(np.array(cube_data)   , axis=0))
	cube_data_clp_med = np.asarray(np.nanmedian(np.array(cube_data), axis=0))
	cube_data_clp_avg = np.asarray(np.nanmean(np.array(cube_data)  , axis=0))
	surfaces_2b_nms   = ['slc_'+str(slc_nmb+1),
						'clp_sum',
						'clp_med',
						'clp_avg']
	surfaces_2b_stt   = [cube_data_sgl_slc,
						cube_data_clp_sum,
						cube_data_clp_med,
						cube_data_clp_avg]
	surfaces_2b_hdr  = ['N',
						'S',
						'M',
						'A']
					
	STT_TBL=['CSL_SNM',
			'CSL_CBW',
			'CSL_FSM',
			'CSL_FAV',
			'CSL_FMD',
			'CSL_FTV',
			'CSL_FSD',
			'CSL_FMX',
			'CSL_FMN',
			'CSL_F25',
			'CSL_F75',
			'CSL_MX1',
			'CSL_SN1',
			'CSL_MX2',
			'CSL_SN2',
			'CSL_LUM',
			'CSL_LLM']
	
	rt = aptbl.Table()
	rt['STT_TBL'] = STT_TBL
	for surface_stat in range(len(surfaces_2b_stt)):
		VAL_TBL=[]
		surface_clpsd   = np.ravel(surfaces_2b_stt[surface_stat])
		fluxtot_sum     = np.nansum(surface_clpsd)
		fluxtot_std     = np.nanstd(surface_clpsd)
		flux_max        = surfaces_2b_stt[surface_stat][x_num/2,y_num/2]
		flux_max_reg    = np.nanmax(surfaces_2b_stt[surface_stat])
		indx_max_reg    = np.where(surfaces_2b_stt[surface_stat] == flux_max_reg)

		cbstt_cbeslnm   = (slc_nmb+1)                                        #Slice Numbr
		cbstt_cbewdth   = (cubewdthv)                                        #Cube Width 
		cbstt_cbeflxt   = (fluxtot_sum)                                      #fluxtot_sum
		cbstt_cbemean   = (np.nanmean(surface_clpsd))                        #fluxtot_avg
		cbstt_cbemedn   = (np.nanmedian(surface_clpsd))                      #fluxtot_med
		cbstt_cbevarc   = (np.nanvar(surface_clpsd))                         #fluxtot_var
		cbstt_cbestdv   = (fluxtot_std)                                      #fluxtot_std
		cbstt_cbemaxm   = (np.nanmax(surface_clpsd))                         #fluxtot_max
		cbstt_cbeminm   = (np.nanmin(surface_clpsd))                         #fluxtot_min
		cbstt_cbept25   = (np.nanpercentile(np.array(surface_clpsd),25))     #fluxtot_p25
		cbstt_cbept75   = (np.nanpercentile(np.array(surface_clpsd),75))     #fluxtot_p75
		cbstt_cbeMAX1   = (flux_max)                                         #MAX_CTR    
		cbstt_cbeMXS1   = (flux_max/fluxtot_std)                             #SNR_MCR    
		cbstt_cbeMAX2   = (flux_max_reg)                                     #MAX_REG    
		cbstt_cbeMXS2   = (flux_max_reg/fluxtot_std)                         #SNR_MRR    
		cbstt_cbeLUMT   = (FluxToLum(fluxtot_sum,z_f2l,frq_r))[0]            #luminosity 
		cbstt_cbeLLMT   = (FluxToLum(fluxtot_sum,z_f2l,frq_r))[1]            #log luminosity 

		VAL_TBL.append(cbstt_cbeslnm)                                        #Slice Numbr 'CSL_SNM'
		VAL_TBL.append(cbstt_cbewdth)                                        #Cube Width  'CSL_CBW'
		VAL_TBL.append(cbstt_cbeflxt)                                        #fluxtot_sum 'CSL_FSM'
		VAL_TBL.append(cbstt_cbemean)                                        #fluxtot_avg 'CSL_FAV'
		VAL_TBL.append(cbstt_cbemedn)                                        #fluxtot_med 'CSL_FMD'
		VAL_TBL.append(cbstt_cbevarc)                                        #fluxtot_var 'CSL_FTV'
		VAL_TBL.append(cbstt_cbestdv)                                        #fluxtot_std 'CSL_FSD'
		VAL_TBL.append(cbstt_cbemaxm)                                        #fluxtot_max 'CSL_FMX'
		VAL_TBL.append(cbstt_cbeminm)                                        #fluxtot_min 'CSL_FMN'
		VAL_TBL.append(cbstt_cbept25)                                        #fluxtot_p25 'CSL_F25'
		VAL_TBL.append(cbstt_cbept75)                                        #fluxtot_p75 'CSL_F75'
		VAL_TBL.append(cbstt_cbeMAX1)                                        #MAX_CTR     'CSL_MX1'
		VAL_TBL.append(cbstt_cbeMXS1)                                        #SNR_MCR     'CSL_SN1'
		VAL_TBL.append(cbstt_cbeMAX2)                                        #MAX_REG     'CSL_MX2'
		VAL_TBL.append(cbstt_cbeMXS2)                                        #SNR_MRR     'CSL_SN2'
		VAL_TBL.append(cbstt_cbeLUMT)                                        #luminosity  'CSL_LUM'
		VAL_TBL.append(cbstt_cbeLLMT)                                        #luminosity  'CSL_LLM'
		rt[str(surfaces_2b_nms[surface_stat])] = VAL_TBL

		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SNM',cbstt_cbeslnm   ,header_comment='Cube Stat Slice Slice Numbr ' + str(surfaces_2b_hdr[surface_stat])) #Slice Numbr
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_CBW',cbstt_cbewdth   ,header_comment='Cube Stat Slice Cube Width '  + str(surfaces_2b_hdr[surface_stat])) #Cube Width 
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FSM',cbstt_cbeflxt   ,header_comment='Cube Stat Slice Flxtot_sum '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_sum
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FAV',cbstt_cbemean   ,header_comment='Cube Stat Slice Flxtot_avg '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_avg
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMD',cbstt_cbemedn   ,header_comment='Cube Stat Slice Flxtot_med '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_med
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FTV',cbstt_cbevarc   ,header_comment='Cube Stat Slice Flxtot_var '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_var
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FSD',cbstt_cbestdv   ,header_comment='Cube Stat Slice Flxtot_std '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_std
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMX',cbstt_cbemaxm   ,header_comment='Cube Stat Slice Flxtot_max '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_max
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMN',cbstt_cbeminm   ,header_comment='Cube Stat Slice Flxtot_min '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_min
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_F25',cbstt_cbept25   ,header_comment='Cube Stat Slice Flxtot_p25 '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_p25
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_F75',cbstt_cbept75   ,header_comment='Cube Stat Slice Flxtot_p75 '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_p75
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_MX1',cbstt_cbeMAX1   ,header_comment='Cube Stat Slice MAX Center '  + str(surfaces_2b_hdr[surface_stat])) #MAX_CTR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SN1',cbstt_cbeMXS1   ,header_comment='Cube Stat Slice SNR Center '  + str(surfaces_2b_hdr[surface_stat])) #SNR_MCR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_MX2',cbstt_cbeMAX2   ,header_comment='Cube Stat Slice MAX Region '  + str(surfaces_2b_hdr[surface_stat])) #MAX_REG    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SN2',cbstt_cbeMXS2   ,header_comment='Cube Stat Slice SNR_Region '  + str(surfaces_2b_hdr[surface_stat])) #SNR_MRR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_LUM',cbstt_cbeLUMT   ,header_comment='Cube Stat Slice Lum '         + str(surfaces_2b_hdr[surface_stat])) #luminosity 
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_LLM',cbstt_cbeLLMT   ,header_comment='Cube Stat Slice Lum [Log] '   + str(surfaces_2b_hdr[surface_stat])) #luminosity 

	rt.write(TABLESTATNAME_2, format=tbl_format_opt, overwrite = True)
	rt.write(TABLESTATNAME_1, format='ascii.fixed_width_two_line', overwrite = True)
	print
	print (colored(Cube2Stt,'yellow'))
	print (colored(Message1,'green'))
	print (colored(Message2,'green'))

def Cube_Stat_2D(Cube2Stt,*args,**kwargs):
	dest_dir  = kwargs.get('dest_dir', None)
	autoaxis  = kwargs.get('autoaxis', False)
	verbose   = kwargs.get('verbose' , False)
	epssave   = kwargs.get('epssave' , False)
	showplot  = kwargs.get('showplot', False)
	clp_fnc   = kwargs.get('clp_fnc' ,'sum')

	z_avg     = kwargs.get('z_avg',Header_Get(Cube2Stt,'STZ_AVG'))
	z_med     = kwargs.get('z_med',Header_Get(Cube2Stt,'STZ_MED'))
	frq_r     = kwargs.get('frq_r',Header_Get(Cube2Stt,'RESTFRQ'))
	z_f2l     = z_med
	cubewdthv = kwargs.get('cubewdthv'      ,1)


	x_ref    = kwargs.get('x_ref',0)
	y_ref    = kwargs.get('y_ref',0)
	ap_size  = kwargs.get('ap_size',0)

	cube_data     = np.asarray(apgtdt(Cube2Stt,memmap=False) )
	freq_num,y_num,x_num = cube_data.shape

	slc_nmb  = kwargs.get('slc_nmb', freq_num/2)

	TABLESTATNAME_1 = (str(Cube2Stt).split('.fits')[0]).split('/')[-1] + '-stt.dat'
	TABLESTATNAME_2 = (str(Cube2Stt).split('.fits')[0]).split('/')[-1] + '-stt.'+ tbl_format_opt
	if dest_dir != None:
		TABLESTATNAME_1 = str(dest_dir)  + '/' + TABLESTATNAME_1
		TABLESTATNAME_2 = str(dest_dir)  + '/' + TABLESTATNAME_2
	elif dest_dir == None:
		TABLESTATNAME_1 = tbl_dir_res    + '/' + TABLESTATNAME_1
		TABLESTATNAME_2 = tbl_dir_res    + '/' + TABLESTATNAME_2
	Message1 = 'Cube Stat on Table: ' + TABLESTATNAME_1 + ' slice number : ' + str(slc_nmb)
	Message2 = 'Cube Stat on Table: ' + TABLESTATNAME_2 + ' slice number : ' + str(slc_nmb)

	cube_data_sgl_slc = cube_data[slc_nmb]
	cube_data_clp_sum = np.asarray(np.nansum(np.array(cube_data)   , axis=0))  #np.nansum(np.array(img_stat)   , axis=0) #
	cube_data_clp_med = np.asarray(np.nanmedian(np.array(cube_data), axis=0))  #np.nansum(np.array(img_stat)   , axis=0) #
	cube_data_clp_avg = np.asarray(np.nanmean(np.array(cube_data)  , axis=0))  #np.nansum(np.array(img_stat)   , axis=0) #
	surfaces_2b_nms   = ['slc_'+str(slc_nmb+1),
						'clp_sum',
						'clp_med',
						'clp_avg']
	surfaces_2b_stt   = [cube_data_sgl_slc,
						cube_data_clp_sum,
						cube_data_clp_med,
						cube_data_clp_avg]
	surfaces_2b_hdr  = ['N',
						'S',
						'M',
						'A']
					
	STT_TBL=['CSL_SNM',
			'CSL_CBW',
			'CSL_FSM',
			'CSL_FAV',
			'CSL_FMD',
			'CSL_FTV',
			'CSL_FSD',
			'CSL_FMX',
			'CSL_FMN',
			'CSL_F25',
			'CSL_F75',
			'CSL_MX1',
			'CSL_SN1',
			'CSL_MX2',
			'CSL_SN2',
			'CSL_LUM',
			'CSL_LLM']
	
	rt = aptbl.Table()
	rt['STT_TBL'] = STT_TBL
	for surface_stat in range(len(surfaces_2b_stt)):
		VAL_TBL=[]
		surface_clpsd   = np.ravel(surfaces_2b_stt[surface_stat])
		fluxtot_sum     = np.nansum(surface_clpsd)
		fluxtot_std     = np.nanstd(surface_clpsd)
		flux_max        = surfaces_2b_stt[surface_stat][x_num/2,y_num/2]
		flux_max_reg    = np.nanmax(surfaces_2b_stt[surface_stat])
		indx_max_reg    = np.where(surfaces_2b_stt[surface_stat] == flux_max_reg)

		cbstt_cbeslnm   = (slc_nmb+1)                                        #Slice Numbr
		cbstt_cbewdth   = (cubewdthv)                                        #Cube Width 
		cbstt_cbeflxt   = (fluxtot_sum)                                      #fluxtot_sum
		cbstt_cbemean   = (np.nanmean(surface_clpsd))                        #fluxtot_avg
		cbstt_cbemedn   = (np.nanmedian(surface_clpsd))                      #fluxtot_med
		cbstt_cbevarc   = (np.nanvar(surface_clpsd))                         #fluxtot_var
		cbstt_cbestdv   = (fluxtot_std)                                      #fluxtot_std
		cbstt_cbemaxm   = (np.nanmax(surface_clpsd))                         #fluxtot_max
		cbstt_cbeminm   = (np.nanmin(surface_clpsd))                         #fluxtot_min
		cbstt_cbept25   = (np.nanpercentile(np.array(surface_clpsd),25))     #fluxtot_p25
		cbstt_cbept75   = (np.nanpercentile(np.array(surface_clpsd),75))     #fluxtot_p75
		cbstt_cbeMAX1   = (flux_max)                                         #MAX_CTR    
		cbstt_cbeMXS1   = (flux_max/fluxtot_std)                             #SNR_MCR    
		cbstt_cbeMAX2   = (flux_max_reg)                                     #MAX_REG    
		cbstt_cbeMXS2   = (flux_max_reg/fluxtot_std)                         #SNR_MRR    
		cbstt_cbeLUMT   = (FluxToLum(fluxtot_sum,z_f2l,frq_r))[0]            #luminosity 
		cbstt_cbeLLMT   = (FluxToLum(fluxtot_sum,z_f2l,frq_r))[1]            #log luminosity 

		VAL_TBL.append(cbstt_cbeslnm)                                        #Slice Numbr 'CSL_SNM'
		VAL_TBL.append(cbstt_cbewdth)                                        #Cube Width  'CSL_CBW'
		VAL_TBL.append(cbstt_cbeflxt)                                        #fluxtot_sum 'CSL_FSM'
		VAL_TBL.append(cbstt_cbemean)                                        #fluxtot_avg 'CSL_FAV'
		VAL_TBL.append(cbstt_cbemedn)                                        #fluxtot_med 'CSL_FMD'
		VAL_TBL.append(cbstt_cbevarc)                                        #fluxtot_var 'CSL_FTV'
		VAL_TBL.append(cbstt_cbestdv)                                        #fluxtot_std 'CSL_FSD'
		VAL_TBL.append(cbstt_cbemaxm)                                        #fluxtot_max 'CSL_FMX'
		VAL_TBL.append(cbstt_cbeminm)                                        #fluxtot_min 'CSL_FMN'
		VAL_TBL.append(cbstt_cbept25)                                        #fluxtot_p25 'CSL_F25'
		VAL_TBL.append(cbstt_cbept75)                                        #fluxtot_p75 'CSL_F75'
		VAL_TBL.append(cbstt_cbeMAX1)                                        #MAX_CTR     'CSL_MX1'
		VAL_TBL.append(cbstt_cbeMXS1)                                        #SNR_MCR     'CSL_SN1'
		VAL_TBL.append(cbstt_cbeMAX2)                                        #MAX_REG     'CSL_MX2'
		VAL_TBL.append(cbstt_cbeMXS2)                                        #SNR_MRR     'CSL_SN2'
		VAL_TBL.append(cbstt_cbeLUMT)                                        #luminosity  'CSL_LUM'
		VAL_TBL.append(cbstt_cbeLLMT)                                        #luminosity  'CSL_LLM'
		rt[str(surfaces_2b_nms[surface_stat])] = VAL_TBL

		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SNM',cbstt_cbeslnm   ,header_comment='Cube Stat Slice Slice Numbr ' + str(surfaces_2b_hdr[surface_stat])) #Slice Numbr
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_CBW',cbstt_cbewdth   ,header_comment='Cube Stat Slice Cube Width '  + str(surfaces_2b_hdr[surface_stat])) #Cube Width 
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FSM',cbstt_cbeflxt   ,header_comment='Cube Stat Slice Flxtot_sum '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_sum
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FAV',cbstt_cbemean   ,header_comment='Cube Stat Slice Flxtot_avg '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_avg
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMD',cbstt_cbemedn   ,header_comment='Cube Stat Slice Flxtot_med '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_med
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FTV',cbstt_cbevarc   ,header_comment='Cube Stat Slice Flxtot_var '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_var
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FSD',cbstt_cbestdv   ,header_comment='Cube Stat Slice Flxtot_std '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_std
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMX',cbstt_cbemaxm   ,header_comment='Cube Stat Slice Flxtot_max '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_max
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_FMN',cbstt_cbeminm   ,header_comment='Cube Stat Slice Flxtot_min '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_min
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_F25',cbstt_cbept25   ,header_comment='Cube Stat Slice Flxtot_p25 '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_p25
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_F75',cbstt_cbept75   ,header_comment='Cube Stat Slice Flxtot_p75 '  + str(surfaces_2b_hdr[surface_stat])) #flxtot_p75
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_MX1',cbstt_cbeMAX1   ,header_comment='Cube Stat Slice MAX Center '  + str(surfaces_2b_hdr[surface_stat])) #MAX_CTR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SN1',cbstt_cbeMXS1   ,header_comment='Cube Stat Slice SNR Center '  + str(surfaces_2b_hdr[surface_stat])) #SNR_MCR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_MX2',cbstt_cbeMAX2   ,header_comment='Cube Stat Slice MAX Region '  + str(surfaces_2b_hdr[surface_stat])) #MAX_REG    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_SN2',cbstt_cbeMXS2   ,header_comment='Cube Stat Slice SNR_Region '  + str(surfaces_2b_hdr[surface_stat])) #SNR_MRR    
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_LUM',cbstt_cbeLUMT   ,header_comment='Cube Stat Slice Lum '         + str(surfaces_2b_hdr[surface_stat])) #luminosity 
		Header_Add(Cube2Stt,'CS'+str(surfaces_2b_hdr[surface_stat])+'_LLM',cbstt_cbeLLMT   ,header_comment='Cube Stat Slice Lum [Log] '   + str(surfaces_2b_hdr[surface_stat])) #luminosity 

	rt.write(TABLESTATNAME_2, format=tbl_format_opt, overwrite = True)
	rt.write(TABLESTATNAME_1, format='ascii.fixed_width_two_line', overwrite = True)
	print
	print (colored(Cube2Stt,'yellow'))
	print (colored(Message1,'green'))
	print (colored(Message2,'green'))
####Fnc_Stk_Stt###