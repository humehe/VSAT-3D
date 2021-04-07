import pandas as pd
from astropy import table as aptbl

from Fnc_Stk_Dir import *

####Fnc_Stk_Tbl####
def split_variable_vars(split_variable):
	if split_variable       == 'RDS' or split_variable == 'RDS_B':
		Splt_Col     = 8
		Splt_Hdr_Cmt = 'Redshift'
		Splt_CNm     = 8
		Splt_Hdr     = 'RDS'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'z'
	elif split_variable     == 'STM' or split_variable == 'STM_B':
		Splt_Col     = 24
		Splt_Hdr_Cmt = 'Stellar Mass [log(M_*/M_sun)]'
		Splt_CNm     = 24
		Splt_Hdr     = 'STM'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'log[$M/M_{\odot}$]'
	elif split_variable     == 'SFR' or split_variable == 'SFR_B':
		Splt_Col     = 22
		Splt_Hdr_Cmt = 'SFR [M_sun/yr]'
		Splt_CNm     = 22
		Splt_Hdr     = 'SFR'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'SFR'
	elif split_variable     == 'LCO' or split_variable == 'LCO_B':
		Splt_Col     = 35
		Splt_Hdr_Cmt = 'CO Lum [K km/s/pc2]'
		Splt_CNm     = 35
		Splt_Hdr     = 'LCO'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'L$_{CO}$'
	elif split_variable     == 'sSF' or split_variable == 'sSF_B':
		Splt_Col     = 26
		Splt_Hdr_Cmt = 'Specific SFR [1/Gyr]'
		Splt_CNm     = 26
		Splt_Hdr     = 'sSF'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'sSFR'
	elif split_variable     == 'MH2' or split_variable == 'MH2_B':
		Splt_Col     = 37
		Splt_Hdr_Cmt = 'H2 mass [log(M_H2/M_sun)]'
		Splt_CNm     = 37
		Splt_Hdr     = 'MH2'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'M$_{H2}$'
	elif split_variable     == 'SFE' or split_variable == 'SFE_B':
		Splt_Col     = 41
		Splt_Hdr_Cmt = 'SFE [1/Gyr]'
		Splt_CNm     = 41
		Splt_Hdr     = 'SFE'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'SFE'
	elif split_variable     == 'LIR' or split_variable == 'LIR_B':
		Splt_Col     = 20
		Splt_Hdr_Cmt = 'LIR [log(L_IR/L_sun)]'
		Splt_CNm     = 20
		Splt_Hdr     = 'LIR'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'log[L$_{IR}$/L${_sun}$]'
	elif split_variable     == 'LFIR' or split_variable == 'LFIR_B':
		Splt_Col     = 11
		Splt_Hdr_Cmt = 'LFIR [log(Lfir/Lo)]'
		Splt_CNm     = 11
		Splt_Hdr     = 'LFIR'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'LFIR [log(Lfir/Lo)]'
	elif split_variable     == 'SDG' or split_variable == 'SDG_B':
		Splt_Col     = 43
		Splt_Hdr_Cmt = 'Surf Dens Gas [log(M_sun/pc2)]'
		Splt_CNm     = 43
		Splt_Hdr     = 'SDG'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = '$\Sigma_{Gas}$'
	elif split_variable     == 'SDS' or split_variable == 'SDS_B':
		Splt_Col     = 45
		Splt_Hdr_Cmt = 'Surf Dens SFR [log(M_sun/yr/kpc2)]'
		Splt_CNm     = 45
		Splt_Hdr     = 'SDS'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = '$\Sigma_{SFR}$'
	elif split_variable     == 'TDT' or split_variable == 'TDT_B':
		Splt_Col     = 47
		Splt_Hdr_Cmt = 'Depletion Time [Gyr]'
		Splt_CNm     = 47
		Splt_Hdr     = 'TDT'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = '$\\tau$'
	elif split_variable       == 'MRP' or split_variable == 'MRP_B':
		Splt_Col     = 50
		Splt_Hdr_Cmt = 'Morphology'	
		Splt_CNm     = 50
		Splt_Hdr     = 'MRP'
		Splt_Hdr_Plt = split_variable
		Splt_Plt_lbl = 'Morphology'
	else:
		print
		print ('split_variable_vars')
		print (colored('Variable '+ str(split_variable) + ' does not exist!','yellow'))
		print
		quit()
	return Splt_Col,Splt_Hdr,Splt_Hdr_Cmt,Splt_CNm,Splt_Hdr_Plt,Splt_Plt_lbl

def Table_Read_org(table_name,format_tbl,*args, **kwargs):
	ftbl = aptbl.Table.read(table_name, format=format_tbl)
	c1   = ftbl['ID']
	c2   = ftbl['fits']
	c3   = ftbl['Source']
	c4   = ftbl['Delta_nu']
	c5   = ftbl['RMS']
	c6   = ftbl['SPW']
	c7   = ftbl['State']
	c8   = ftbl['z']
	c9   = ftbl['RA']
	c10  = ftbl['Dec']
	c11  = ftbl['log(Lfir/Lo)']
	c12  = ftbl['D_log(Lfir/Lo)']
	c13  = ftbl['nu_obs']
	c14  = ftbl['V_obs']
	c15  = ftbl['Vobs_err']
	c16  = ftbl['DV']
	c17  = ftbl['DV_err']
	c18  = ftbl['Maxis_2']
	c19  = ftbl['Maxis_err_2']
	c20  = ftbl['maxis_2a']
	c21  = ftbl['maxis_err_2a']
	c22  = ftbl['angle']
	c23  = ftbl['angle_err']
	c24  = ftbl['COFlux']
	c25  = ftbl['err_COflux']
	c26  = ftbl['D_L']
	c27  = ftbl['D_A']
	c28  = ftbl['logMstar']
	c29  = ftbl['er_logMstar']
	c30  = ftbl['Tam_Fit']
	c31  = ftbl['Bean_size']
	c32  = ftbl['Morfology_Coments']
	c33  = ftbl['R_petro_r']
	c34  = ftbl['g']
	c35  = ftbl['r']
	c36  = ftbl['mu_i']
	c37  = ftbl['H_i']
	return(ftbl,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)

def Table_Read(table_name,format_tbl,*args, **kwargs):
	ftbl = aptbl.Table.read(table_name, format=format_tbl)
	c1   = ftbl['ID']
	c2   = ftbl['fits']
	c3   = ftbl['Source']
	c4   = ftbl['Delta_nu']
	c5   = ftbl['RMS']
	c6   = ftbl['SPW']
	c7   = ftbl['State']
	c8   = ftbl['z_1']
	c9   = ftbl['RA']
	c10  = ftbl['Dec']
	c11  = ftbl['log(Lfir/Lo)']
	c12  = ftbl['D_log(Lfir/Lo)']
	c13  = ftbl['nu_obs']
	c14  = ftbl['V_obs']
	c15  = ftbl['GAMA_ID']
	c16  = ftbl['SOURCE']
	c17  = ftbl['RAJ2000']
	c18  = ftbl['DECJ2000']
	c19  = ftbl['z_2']
	c20  = ftbl['log[L_IR/L_sun]']
	c21  = ftbl['c7_err']
	c22  = ftbl['SFR']
	c23  = ftbl['c8_err']
	c24  = ftbl['log[M_S/M_sun]']
	c25  = ftbl['c9_err']
	c26  = ftbl['sSFR']
	c27  = ftbl['c10_err']
	c28  = ftbl['nu_ob_1']
	c29  = ftbl['nu_ob_2']
	c30  = ftbl['c11_err']
	c31  = ftbl['v_fwhm']
	c32  = ftbl['c12_err']
	c33  = ftbl['S_COXDeltaV']
	c34  = ftbl['c13_err']
	c35  = ftbl['L_CO']
	c36  = ftbl['c14_err']
	c37  = ftbl['log[M_H2/M_sun]']
	c38  = ftbl['c15_err']
	c39  = ftbl['R_FWHM']
	c40  = ftbl['c16_err']
	c41  = ftbl['SFE']
	c42  = ftbl['c17_err']
	c43  = ftbl['SIGMA_gas']
	c44  = ftbl['c18_err']
	c45  = ftbl['SIGMA_SFR']
	c46  = ftbl['c19_err']
	c47  = ftbl['Tau_gas']
	c48  = ftbl['c20_err']
	if 'MRP' in table_name:
		c49  = ftbl['M']
		c50  = ftbl['M_N']
		return(ftbl,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
			c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
			c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,
			c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,
			c41,c42,c43,c44,c45,c46,c47,c48,c49,c50)
	else:
		return(ftbl,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
			c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
			c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,
			c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,
			c41,c42,c43,c44,c45,c46,c47,c48)

def Table_Read_Lit(table_name,format_tbl,*args, **kwargs):
	ftbl = aptbl.Table.read(table_name, format=format_tbl)
	if 'SDS' in table_name:
		c1  = ftbl['Name']
		c2  = ftbl['SSFR']
		c4  = ftbl['R1']
		c5  = ftbl['12CO-13CO']
		c6  = ftbl['12CO-13CO_e']
		c7  = ftbl['R2']
		c8  = ftbl['12CO-13CO-2']
		c9  = ftbl['12CO-13CO-2_e']
		c10 = ftbl['R3']
		c11 = ftbl['12CO/C18O']
		c12 = ftbl['12CO/C18O_e']
		c13 = ftbl['R4']
		return(ftbl,c1,c2,c5,c6,c4,c7,c8,c9,c10,
			c11,c12,c13)
	else:
		pass

def Table_Ipt_Cat_Stats(Cat_Ipt_Tbl_Sts,hdr_sts):
	Splt_Vars = split_variable_vars(hdr_sts)
	Tbl_Splt_Col = Splt_Vars[0]
	Tbl_Splt_Hdr = Splt_Vars[1]
	Tbl_Splt_Hdr_Cmt = Splt_Vars[2]

	z_sample_avg     = np.mean(Cat_Ipt_Tbl_Sts[8])
	z_sample_med     = np.median(Cat_Ipt_Tbl_Sts[8])
	z_sample_1sl     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 15.9)
	z_sample_1sh     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 84.1)
	z_sample_2sl     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 2.30)
	z_sample_2sh     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 97.7)
	z_sample_3sl     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 0.20)
	z_sample_3sh     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 99.8)
	z_sample_p25     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 25.0)
	z_sample_p75     = np.nanpercentile(Cat_Ipt_Tbl_Sts[8], 75.0)

	Splt_sample_avg  = np.mean(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col])
	Splt_sample_med  = np.median(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col])
	Splt_sample_1sl  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 15.9)
	Splt_sample_1sh  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 84.1)
	Splt_sample_2sl  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 2.30)
	Splt_sample_2sh  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 97.7)
	Splt_sample_3sl  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 0.20)
	Splt_sample_3sh  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 99.8)
	Splt_sample_p25  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 25.0)
	Splt_sample_p75  = np.nanpercentile(Cat_Ipt_Tbl_Sts[Tbl_Splt_Col], 75.0)
	print
	print ('Redshift (avg): ',z_sample_avg)
	print ('Redshift (med): ',z_sample_med)
	print ('Redshift (avg): ',Splt_sample_avg)
	print ('Redshift (med): ',Splt_sample_med)
	print ('subcube_width : ',subcube_width)
	var_sts = [
				['STZ_AVG','STZ_MED',
				'STZ_1SL','STZ_1SH',
				'STZ_2SL','STZ_2SH',
				'STZ_3SL','STZ_3SH',
				'STZ_P25','STZ_P75',
				'STS_AVG','STS_MED',
				'STS_1SL','STS_1SH',
				'STS_2SL','STS_2SH',
				'STS_3SL','STS_3SH',
				'STS_P25','STS_P75'],
				[z_sample_avg,z_sample_med,
				z_sample_1sl,z_sample_1sh,
				z_sample_2sl,z_sample_2sh,
				z_sample_3sl,z_sample_3sh,
				z_sample_p25,z_sample_p75,
				Splt_sample_avg,Splt_sample_med,
				Splt_sample_1sl,Splt_sample_1sh,
				Splt_sample_2sl,Splt_sample_2sh,
				Splt_sample_3sl,Splt_sample_3sh,
				Splt_sample_p25,Splt_sample_p75],
				Tbl_Splt_Hdr,
				Tbl_Splt_Hdr_Cmt					
				]
	return var_sts

def Read_MCMC_Table_Sample_Stat(tbl_mcmc_ipt,format_tbl,tbl_typ,line,*args, **kwargs):
	ftbl = aptbl.Table.read(tbl_mcmc_ipt, format=format_tbl)
	if tbl_typ == 'Z' and line == 1:
		c1  = ftbl.columns[1] #1-N1_TOT 		9-N2_TOT
		c2  = ftbl.columns[2] #2-Z1_MED 		10-Z2_MED
		c3  = ftbl.columns[3] #3-Z1_MED_E1SGL 	11-Z2_MED_E1SGL
		c4  = ftbl.columns[4] #4-Z1_MED_E1SGH 	12-Z2_MED_E1SGH
		return ftbl,c1,c2,c3,c4
	elif tbl_typ == 'Z' and line == 2:
		c1  = ftbl.columns[9]  #9-N2_TOT
		c2  = ftbl.columns[10] #10-Z2_MED
		c3  = ftbl.columns[11] #11-Z2_MED_E1SGL
		c4  = ftbl.columns[12] #12-Z2_MED_E1SGH
		return ftbl,c1,c2,c3,c4
	elif tbl_typ == 'Flx' and line == 1:
		c1  = ftbl.columns[1]	#1-N1_TOT
		c2  = ftbl.columns[2]	#2-S1_AVG
		c3  = ftbl.columns[3]	#3-S1_AVG_E
		c4  = ftbl.columns[4]	#4-S1_MED
		c5  = ftbl.columns[5]	#5-S1_MED_E
		return ftbl,c1,c2,c3,c4,c5
	elif tbl_typ == 'Flx' and line == 2:
		c1  = ftbl.columns[6]	#6-N1_TOT
		c2  = ftbl.columns[7]	#7-S1_AVG
		c3  = ftbl.columns[8]	#8-S1_AVG_E
		c4  = ftbl.columns[9]	#9-S1_MED
		c5  = ftbl.columns[10]	#10-S1_MED_E
		return ftbl,c1,c2,c3,c4,c5
	elif tbl_typ == 'Lum' and line == 1:
		c1  = ftbl.columns[1]	#N1_TOT
		c2  = ftbl.columns[2]	#L1_AVG
		c3  = ftbl.columns[3]	#L1_AVG_E1SGL
		c4  = ftbl.columns[4]	#L1_AVG_E1SGH
		c10 = ftbl.columns[9]	#L1_MED
		c11 = ftbl.columns[10]	#L1_MED_E1SGL
		c12 = ftbl.columns[11]	#L1_MED_E1SGH
		return ftbl,c1,c2,c3,c4,c10,c11,c12
	elif tbl_typ == 'Lum' and line == 2:
		c1  = ftbl.columns[17]	#N1_TOT
		c2  = ftbl.columns[18]	#L2_AVG
		c3  = ftbl.columns[19]	#L2_AVG_E1SGL
		c4  = ftbl.columns[20]	#L2_AVG_E1SGH
		c10 = ftbl.columns[25]	#L2_MED
		c11 = ftbl.columns[26]	#L2_MED_E1SGL
		c12 = ftbl.columns[27]	#L2_MED_E1SGH
		return ftbl,c1,c2,c3,c4,c10,c11,c12
	elif tbl_typ == 'Var' and line ==1:
		c1  = ftbl.columns[1]
		c2  = ftbl.columns[2]
		c3  = ftbl.columns[3]
		c4  = ftbl.columns[4]
		return ftbl,c1,c2,c3,c4
	elif tbl_typ == 'Var' and line ==2:
		c1  = ftbl.columns[1]
		c2  = ftbl.columns[2]
		c3  = ftbl.columns[3]
		c4  = ftbl.columns[4]
		return ftbl,c1,c2,c3,c4

	else:
		pass
	
def Read_MCMC_Table_MCMC_PLT(tbl_mcmc_plt_ipt,format_tbl,*args,**kwargs):
	ftbl = aptbl.Table.read(tbl_mcmc_plt_ipt, format=format_tbl)
	Z1_AVG_E1_MC_OPT = ftbl['Z1_AVG_E1_MC_OPT'] 
	Z2_AVG_E1_MC_OPT = ftbl['Z2_AVG_E1_MC_OPT'] 
	Z1_MED_E1_MC_OPT = ftbl['Z1_MED_E1_MC_OPT'] 
	Z2_MED_E1_MC_OPT = ftbl['Z2_MED_E1_MC_OPT'] 
	V1_AVG_E1_MC_OPT = ftbl['V1_AVG_E1_MC_OPT'] 
	V2_AVG_E1_MC_OPT = ftbl['V2_AVG_E1_MC_OPT'] 
	V1_MED_E1_MC_OPT = ftbl['V1_MED_E1_MC_OPT'] 
	V2_MED_E1_MC_OPT = ftbl['V2_MED_E1_MC_OPT'] 
	S1_AVG_E1_MC_OPT = ftbl['S1_AVG_E1_MC_OPT'] 
	S2_AVG_E1_MC_OPT = ftbl['S2_AVG_E1_MC_OPT'] 
	S1_MED_E1_MC_OPT = ftbl['S1_MED_E1_MC_OPT'] 
	S2_MED_E1_MC_OPT = ftbl['S2_MED_E1_MC_OPT'] 
	L1_AVG_E1_MC_1_OPT = ftbl['L1_AVG_E1_MC_1_OPT'] 
	L1_MED_E1_MC_1_OPT = ftbl['L1_MED_E1_MC_1_OPT'] 
	L1_AVG_E1_MC_2_OPT = ftbl['L1_AVG_E1_MC_2_OPT'] 
	L1_MED_E1_MC_2_OPT = ftbl['L1_MED_E1_MC_2_OPT'] 
	L2_AVG_E1_MC_1_OPT = ftbl['L2_AVG_E1_MC_1_OPT'] 
	L2_MED_E1_MC_1_OPT = ftbl['L2_MED_E1_MC_1_OPT'] 
	L2_AVG_E1_MC_2_OPT = ftbl['L2_AVG_E1_MC_2_OPT'] 
	L2_MED_E1_MC_2_OPT = ftbl['L2_MED_E1_MC_2_OPT'] 
	S12_AVG_E1_MC_OPT = ftbl['S12_AVG_E1_MC_OPT'] 
	S12_MED_E1_MC_OPT = ftbl['S12_MED_E1_MC_OPT'] 
	L12_AVG_E1_MC_1_OPT = ftbl['L12_AVG_E1_MC_1_OPT'] 
	L12_MED_E1_MC_1_OPT = ftbl['L12_MED_E1_MC_1_OPT'] 
	L12_AVG_E1_MC_2_OPT = ftbl['L12_AVG_E1_MC_2_OPT'] 
	L12_MED_E1_MC_2_OPT = ftbl['L12_MED_E1_MC_2_OPT'] 

	Z1_AVG_E2_MC_OPT = ftbl['Z1_AVG_E2_MC_OPT'] 
	Z2_AVG_E2_MC_OPT = ftbl['Z2_AVG_E2_MC_OPT'] 
	Z1_MED_E2_MC_OPT = ftbl['Z1_MED_E2_MC_OPT'] 
	Z2_MED_E2_MC_OPT = ftbl['Z2_MED_E2_MC_OPT'] 
	V1_AVG_E2_MC_OPT = ftbl['V1_AVG_E2_MC_OPT'] 
	V2_AVG_E2_MC_OPT = ftbl['V2_AVG_E2_MC_OPT'] 
	V1_MED_E2_MC_OPT = ftbl['V1_MED_E2_MC_OPT'] 
	V2_MED_E2_MC_OPT = ftbl['V2_MED_E2_MC_OPT'] 
	S1_AVG_E2_MC_OPT = ftbl['S1_AVG_E2_MC_OPT'] 
	S2_AVG_E2_MC_OPT = ftbl['S2_AVG_E2_MC_OPT'] 
	S1_MED_E2_MC_OPT = ftbl['S1_MED_E2_MC_OPT'] 
	S2_MED_E2_MC_OPT = ftbl['S2_MED_E2_MC_OPT'] 
	L1_AVG_E2_MC_1_OPT = ftbl['L1_AVG_E2_MC_1_OPT'] 
	L1_MED_E2_MC_1_OPT = ftbl['L1_MED_E2_MC_1_OPT'] 
	L1_AVG_E2_MC_2_OPT = ftbl['L1_AVG_E2_MC_2_OPT'] 
	L1_MED_E2_MC_2_OPT = ftbl['L1_MED_E2_MC_2_OPT'] 
	L2_AVG_E2_MC_1_OPT = ftbl['L2_AVG_E2_MC_1_OPT'] 
	L2_MED_E2_MC_1_OPT = ftbl['L2_MED_E2_MC_1_OPT'] 
	L2_AVG_E2_MC_2_OPT = ftbl['L2_AVG_E2_MC_2_OPT'] 
	L2_MED_E2_MC_2_OPT = ftbl['L2_MED_E2_MC_2_OPT'] 
	S12_AVG_E2_MC_OPT = ftbl['S12_AVG_E2_MC_OPT'] 
	S12_MED_E2_MC_OPT = ftbl['S12_MED_E2_MC_OPT'] 
	L12_AVG_E2_MC_1_OPT = ftbl['L12_AVG_E2_MC_1_OPT'] 
	L12_MED_E2_MC_1_OPT = ftbl['L12_MED_E2_MC_1_OPT'] 
	L12_AVG_E2_MC_2_OPT = ftbl['L12_AVG_E2_MC_2_OPT'] 
	L12_MED_E2_MC_2_OPT	 = ftbl['L12_MED_E2_MC_2_OPT'] 
	opt_tbl = [ftbl,
				Z1_AVG_E1_MC_OPT,Z2_AVG_E1_MC_OPT,Z1_MED_E1_MC_OPT,Z2_MED_E1_MC_OPT, #1
				V1_AVG_E1_MC_OPT,V2_AVG_E1_MC_OPT,V1_MED_E1_MC_OPT,V2_MED_E1_MC_OPT, #5
				S1_AVG_E1_MC_OPT,S2_AVG_E1_MC_OPT,S1_MED_E1_MC_OPT,S2_MED_E1_MC_OPT, #9
				L1_AVG_E1_MC_1_OPT,L1_MED_E1_MC_1_OPT,L1_AVG_E1_MC_2_OPT,L1_MED_E1_MC_2_OPT, #13
				L2_AVG_E1_MC_1_OPT,L2_MED_E1_MC_1_OPT,L2_AVG_E1_MC_2_OPT,L2_MED_E1_MC_2_OPT, #17
				S12_AVG_E1_MC_OPT,S12_MED_E1_MC_OPT, #21
				L12_AVG_E1_MC_1_OPT,L12_MED_E1_MC_1_OPT,L12_AVG_E1_MC_2_OPT,L12_MED_E1_MC_2_OPT, #23
				Z1_AVG_E2_MC_OPT,Z2_AVG_E2_MC_OPT,Z1_MED_E2_MC_OPT,Z2_MED_E2_MC_OPT, #27
				V1_AVG_E2_MC_OPT,V2_AVG_E2_MC_OPT,V1_MED_E2_MC_OPT,V2_MED_E2_MC_OPT, #31
				S1_AVG_E2_MC_OPT,S2_AVG_E2_MC_OPT,S1_MED_E2_MC_OPT,S2_MED_E2_MC_OPT, #35
				L1_AVG_E2_MC_1_OPT,L1_MED_E2_MC_1_OPT,L1_AVG_E2_MC_2_OPT,L1_MED_E2_MC_2_OPT, #39
				L2_AVG_E2_MC_1_OPT,L2_MED_E2_MC_1_OPT,L2_AVG_E2_MC_2_OPT,L2_MED_E2_MC_2_OPT, #43
				S12_AVG_E2_MC_OPT,S12_MED_E2_MC_OPT, #47
				L12_AVG_E2_MC_1_OPT,L12_MED_E2_MC_1_OPT,L12_AVG_E2_MC_2_OPT,L12_MED_E2_MC_2_OPT #49
			]
	return opt_tbl
####Fnc_Stk_Tbl####