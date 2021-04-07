import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib import colors as mcolors
import matplotlib.ticker as mticker
import scipy.integrate as integrate

from matplotlib.backends.backend_pdf import PdfPages

from Fnc_Stk_Dir import *
from Fnc_Stk_Fts import *
from Fnc_Stk_Spc import *
from Fnc_Stk_Mth import *
from Fnc_Stk_Tbl import *

####Fnc_Stk_Plt###
def align_yaxis(ax1, v1, ax2, v2):
	"""adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
	_, y1 = ax1.transData.transform((0, v1))
	_, y2 = ax2.transData.transform((0, v2))
	inv = ax2.transData.inverted()
	_, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
	miny, maxy = ax2.get_ylim()
	ax2.set_ylim(miny+dy, maxy+dy)

def align_xaxis(ax1, v1, ax2, v2, y2min, y2max):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1."""

    """where y2max is the maximum value in your secondary plot. I haven't
     had a problem with minimum values being cut, so haven't set this. This
     approach doesn't necessarily make for axis limits at nice near units,
     but does optimist plot space"""

    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    scale = 1
    while scale*(maxy+dy) < y2max:
        scale += 0.05

class ScaledLocator(mpl.ticker.MaxNLocator):
    """
    Locates regular intervals along an axis scaled by *dx* and shifted by
    *x0*. For example, this would locate minutes on an axis plotted in seconds
    if dx=60.  This differs from MultipleLocator in that an approriate interval
    of dx units will be chosen similar to the default MaxNLocator.
    """
    def __init__(self, dx=1.0, x0=0.0):
        self.dx = dx
        self.x0 = x0
        mpl.ticker.MaxNLocator.__init__(self, nbins=9, steps=[1, 2, 5, 10])

    def rescale(self, x):
        return x / self.dx + self.x0
    def inv_rescale(self, x):
        return  (x - self.x0) * self.dx

    def __call__(self): 
        vmin, vmax = self.axis.get_view_interval()
        vmin, vmax = self.rescale(vmin), self.rescale(vmax)
        vmin, vmax = mpl.transforms.nonsingular(vmin, vmax, expander = 0.05)
        locs = self.bin_boundaries(vmin, vmax)
        locs = self.inv_rescale(locs)
        prune = self._prune
        if prune=='lower':
            locs = locs[1:]
        elif prune=='upper':
            locs = locs[:-1]
        elif prune=='both':
            locs = locs[1:-1]
        return self.raise_if_exceeds(locs)

class ScaledFormatter(mpl.ticker.OldScalarFormatter):
    """Formats tick labels scaled by *dx* and shifted by *x0*."""
    def __init__(self, dx=1.0, x0=0.0, **kwargs):
        self.dx, self.x0 = dx, x0

    def rescale(self, x):
        return x / self.dx + self.x0

    def __call__(self, x, pos=None):
        xmin, xmax = self.axis.get_view_interval()
        xmin, xmax = self.rescale(xmin), self.rescale(xmax)
        d = abs(xmax - xmin)
        x = self.rescale(x)
        s = self.pprint_val(x, d)
        return s

def fmt(x, pos):
    a, b = '{:.1f}e-3'.format(x/1e-3).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def Plot_Cube_2D(Cube2bplot_2D,*args,**kwargs):
    dest_dir_plt = kwargs.get('dest_dir_plt',None)
    dest_dir_clp = kwargs.get('dest_dir_clp',None)
    autoaxis     = kwargs.get('autoaxis',False)
    verbose      = kwargs.get('verbose' , False)
    epssave      = kwargs.get('epssave' , False)
    showplot     = kwargs.get('showplot', False) 
    slc_nmb      = kwargs.get('slc_nmb' , None) 
    clp_fnc      = kwargs.get('clp_fnc' , 'sum')

    redshift     = kwargs.get('redshift' ,'1')
    rst_frq      = kwargs.get('rst_frq'  ,'1')

    x_ref        = kwargs.get('x_ref',0)
    y_ref        = kwargs.get('y_ref',0)
    ap_size      = kwargs.get('ap_size',0)

    z_avg        = kwargs.get('z_avg',Header_Get(Cube2bplot_2D,'STZ_AVG'))
    z_med        = kwargs.get('z_med',Header_Get(Cube2bplot_2D,'STZ_MED'))
    frq_r        = kwargs.get('frq_r',Header_Get(Cube2bplot_2D,'RESTFRQ'))
    z_f2l        = z_med

    prefix       = kwargs.get('prefix','')

    dest_dir_plt = kwargs.get('dest_dir_plt',None)

    Cube_Info    = Cube_Header_Get(Cube2bplot_2D,frq_r* u.Hz)
    FRQ_AXS      = Cube_Info[16].value
    VEL_AXS      = Cube_Info[17].value

    flx_scl      = kwargs.get('flx_scl',1e-06)

    if slc_nmb != None:
        pass
    elif slc_nmb == None:
        slc_nmb = 0

    if flx_scl == 1e-06:
        scl_prfx = '$\mu$'
    elif flx_scl == 1e-03:
        scl_prfx = 'm'
    else :
        scl_prfx = ''

    if dest_dir_plt != None:
        PLOTFILENAME_2DS = str(dest_dir_plt)  + '/' + prefix + (str(Cube2bplot_2D).split('.fits')[0]).split('/')[-1] + '-2DS.pdf'
        PLOTFILENAME_CLP = str(dest_dir_plt)  + '/' + prefix + (str(Cube2bplot_2D).split('.fits')[0]).split('/')[-1] + '-2DC-'+ clp_fnc +'.pdf'
    elif dest_dir_plt == None:
        PLOTFILENAME_2DS = stm_dir_plt        + '/' + prefix + (str(Cube2bplot_2D).split('.fits')[0]).split('/')[-1] + '-2DS.pdf'
        PLOTFILENAME_CLP = stm_dir_plt        + '/' + prefix + (str(Cube2bplot_2D).split('.fits')[0]).split('/')[-1] + '-2DC-'+ clp_fnc +'.pdf'

    if dest_dir_clp != None:
        Cube2bclp_2D_opt = dest_dir_clp + (Cube2bplot_2D.split('.fits',1)[0]).rsplit('/',1)[-1] + '-2DC-'+clp_fnc+'.fits'
    elif dest_dir_clp == None:
        Cube2bclp_2D_opt = stp_dir_res + (Cube2bplot_2D.split('.fits',1)[0]).rsplit('/',1)[-1] + '-2DC-'+clp_fnc+'.fits'


    data_2b_plt = np.asarray(apgtdt(Cube2bplot_2D,memmap=False) )

    slice_fwhm = (Header_Get(Cube2bplot_2D,'FTS_FWH'))
    slice_cwdt = (Header_Get(Cube2bplot_2D,'STT_VEL')) 
    slice_nmbr = (Header_Get(Cube2bplot_2D,'MAX_SNS')) 

    slice_wdnb = int(np.ceil(slice_fwhm / slice_cwdt))
    slice_nblw = int(slice_nmbr-int(np.ceil(slice_fwhm / slice_cwdt)))
    slice_nbhg = int(slice_nmbr+int(np.ceil(slice_fwhm / slice_cwdt)))


    if (int(slice_nmbr) == slice_nblw) & (int(slice_nmbr)==slice_nbhg):
        data_2b_plt_clp = data_2b_plt[int(slice_nmbr)]
        data_2d_clp     = data_2b_plt_clp
        Message1        = 'Creating collapsed datacube ('+str(clp_fnc)+')'
        Message2        = 'FWHM      : ' + str(slice_fwhm)
        Message3        = 'Channels  : ' + str(slice_nblw)+'-'+str(slice_nbhg)
        Message4        = 'Fits File : ' + Cube2bclp_2D_opt
        print
        print (colored(Message1,'yellow'))
        print (colored(Message2,'yellow'))
        print
    else:
        data_2b_plt_clp = data_2b_plt[slice_nblw:slice_nbhg]
        if clp_fnc == 'sum':
            data_2d_clp = np.asarray(np.nansum(np.array(data_2b_plt_clp)   , axis=0))
            doublecollapsed = np.ravel(data_2d_clp)
            data_2d_clp =  data_2d_clp
        elif clp_fnc == 'med':
            data_2d_clp = np.asarray(np.nanmedian(np.array(data_2b_plt_clp), axis=0))
            data_2d_clp =  data_2d_clp
            doublecollapsed = np.ravel(data_2d_clp)
            data_2d_clp =  data_2d_clp
        elif clp_fnc == 'avg':
            data_2d_clp = np.asarray(np.nanmean(np.array(data_2b_plt_clp)  , axis=0))
            data_2d_clp =  data_2d_clp
            doublecollapsed = np.ravel(data_2d_clp)
            data_2d_clp =  data_2d_clp

    Wrt_FITS_File(data_2d_clp,Cube2bclp_2D_opt)
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_AVG')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_MED')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_1SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_1SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_2SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_2SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_3SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_3SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_P25')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STZ_P75')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_AVG')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_MED')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_1SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_1SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_2SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_2SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_3SL')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_3SH')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_P25')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'STS_P75')
    Header_Copy(Cube2bclp_2D_opt,Cube2bplot_2D,'RESTFRQ')

    F_plt,X_plt,Y_plt = data_2b_plt.shape
    X_clp,Y_clp       = data_2d_clp.shape

    nx_f2DG, ny_f2DG = X_clp,Y_clp
    nx,ny            = nx_f2DG,ny_f2DG

    X0_f2DG     = kwargs.get('X0_f2DG',nx_f2DG/2)
    Y0_f2DG     = kwargs.get('Y0_f2DG',ny_f2DG/2)
    MAX_CTR     = data_2d_clp[X_plt/2,Y_plt/2]

    ##########################################2DS###################################

    fxsize=9
    fysize=8
    f = plt.figure(num=None, figsize=(fxsize, fysize), dpi=180, facecolor='w',
        edgecolor='k')
    plt.subplots_adjust(
        left    = (25/25.4)/fxsize,       
        bottom  = (16/25.4)/fysize,       
        right   = 1 - (34/25.4)/fxsize,   
        top     = 1 - (15/25.4)/fysize)   
    plt.subplots_adjust(hspace=0)


    gs0 = gridspec.GridSpec(1, 1)
    

    gs11 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
        
    ax110 = plt.Subplot(f, gs11[0,0])
    f.add_subplot(ax110)

    ax110.set_rasterization_zorder(1)
    plt.autoscale(enable=True, axis='y', tight=False)
    ax110.xaxis.set_tick_params(labelsize=16)
    ax110.yaxis.set_tick_params(labelsize=16)
    ax110.set_title('Slice number: '+ str(slc_nmb+1)  +  '-' + str(round(VEL_AXS[slc_nmb],0)) + ' km/s'+
                    "\n" +  'X$_{\mathrm{c}}$,Y$_{\mathrm{c}}$: ' + str(+x_ref) + ','+str(+y_ref) + ' $\\varnothing$: '+ str(ap_size) + '"',
                    family='serif',fontsize=16)
    xticklabels = ax110.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    yticklabels = ax110.get_yticklabels()
    plt.setp(yticklabels, visible=True)

    plt.tick_params(which='both', width=1.0)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=5)
    ax110.minorticks_on()

    plt.xlabel('X',fontsize=16,family = 'serif')
    plt.ylabel('Y',fontsize=16,family = 'serif')
    if ('_ms.' in Cube2bplot_2D) or ('dta_in.' in Cube2bplot_2D) or ('dta_ot.' in Cube2bplot_2D):
        tick_color = 'white'
    elif ('msk_in.' in Cube2bplot_2D) or ('crc.' in Cube2bplot_2D) or ('msk_ot.' in Cube2bplot_2D):
        tick_color = 'black'
    else:
        tick_color = 'white'

    ax110.xaxis.set_tick_params(which='both',labelsize=20,direction='in',color=tick_color,bottom=True,top=True,left=True,right=True)
    ax110.yaxis.set_tick_params(which='both',labelsize=20,direction='in',color=tick_color,bottom=True,top=True,left=True,right=True)

    plt.imshow(data_2b_plt[slc_nmb]/flx_scl, origin='lower',cmap='viridis')
    divider = make_axes_locatable(ax110)
    cax  = divider.append_axes("right", size="5%", pad=0.05)    
    cbar = plt.colorbar(cax=cax,pad=0.5)
    cbar.set_label('S ['+scl_prfx+'Jy]', rotation=270,family = 'serif',size=16,labelpad=15)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.set_yticklabels(["{:.0f}".format(i) for i in cbar.get_ticks()]) 

    ax110.xaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-X0_f2DG+x_ref))
    ax110.yaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-Y0_f2DG+y_ref))

    plt.scatter(X0_f2DG, Y0_f2DG, s=25, c='white', marker='x')
    
    plt.savefig(PLOTFILENAME_2DS)

    Message3       = 'Generated plot file for datacubes slice ('+str(slc_nmb+1)+')'
    Message4       = PLOTFILENAME_2DS
    print
    print (colored(Message3,'cyan'))
    print (colored(Message4,'cyan'))
    print
    ##########################################2DS###################################

    ##########################################CLP###################################

    f = plt.figure(num=None, figsize=(fxsize, fysize), dpi=180, facecolor='w',
        edgecolor='k')
    plt.subplots_adjust(
        left    = (25/25.4)/fxsize,       
        bottom  = (16/25.4)/fysize,       
        right   = 1 - (34/25.4)/fxsize,   
        top     = 1 - (15/25.4)/fysize)   
    plt.subplots_adjust(hspace=0)


    gs0 = gridspec.GridSpec(1, 1)
    

    gs11 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
        
    ax110 = plt.Subplot(f, gs11[0,0])
    f.add_subplot(ax110)

    ax110.set_rasterization_zorder(1)
    plt.autoscale(enable=True, axis='y', tight=False)
    ax110.xaxis.set_tick_params(labelsize=16)
    ax110.yaxis.set_tick_params(labelsize=16)
    ax110.set_title('2D Collapse [' +str(slice_nblw+1) + ':'+ str(slice_nbhg+1)+ ']: ' + str(clp_fnc.upper()) +
                    "\n" +  'X$_{\mathrm{c}}$,Y$_{\mathrm{c}}$: ' + str(+x_ref) + ','+str(+y_ref) + ' $\\varnothing$: '+ str(ap_size) + '"',
                    family='serif',fontsize=16)

    xticklabels = ax110.get_xticklabels()
    plt.setp(xticklabels, visible=True)
    yticklabels = ax110.get_yticklabels()
    plt.setp(yticklabels, visible=True)

    plt.tick_params(which='both', width=1.0)
    plt.tick_params(which='major', length=10)
    plt.tick_params(which='minor', length=5)
    ax110.minorticks_on()

    ax110.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color=tick_color,bottom=True,top=True,left=True,right=True)
    ax110.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color=tick_color,bottom=True,top=True,left=True,right=True)

    plt.xlabel('X',fontsize=16,family = 'serif')
    plt.ylabel('Y',fontsize=16,family = 'serif')
    
    plt.imshow(data_2d_clp/flx_scl, origin='lower',cmap='viridis')
    divider = make_axes_locatable(ax110)
    cax  = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    cbar.set_label('S ['+scl_prfx+'Jy]', rotation=270,family = 'serif',size=16,labelpad=15)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.set_yticklabels(["{:.0f}".format(i) for i in cbar.get_ticks()])

    ax110.xaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-X0_f2DG+x_ref))
    ax110.yaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-Y0_f2DG+y_ref))

    plt.scatter(X0_f2DG, Y0_f2DG, s=25, c='white', marker='x')
    plt.savefig(PLOTFILENAME_CLP)
    Message5       = 'Generated plot file for collapsed datacube ('+str(clp_fnc)+')'
    Message6       = PLOTFILENAME_CLP
    print
    print (colored(Message5,'cyan'))
    print (colored(Message6,'cyan'))
    print

    ##########################################CLP###################################
    if verbose == True:
        print
        print (colored('Generated Fits: ' + str(Cube2bclp_2D_opt) + ' Dim: ' + str(F_plt) + ' X ' + str(X_plt) + ' X ' + str(Y_plt) ,'yellow'))
        print (colored('Generated Plot: ' + str(PLOTFILENAME_2DS) + ' Dim: ' + str(F_plt) + ' X ' + str(X_plt) + ' X ' + str(Y_plt) ,'cyan'))
        print (colored('Generated Plot: ' + str(PLOTFILENAME_CLP) + ' Dim: ' + str(X_clp) + ' X ' + str(Y_clp) ,'cyan'))
    elif verbose ==False:
        pass
    plt.close('all')

def Plot_Cube_3D(Cube2bplot_3D,*args,**kwargs):
    dest_dir  = kwargs.get('dest_dir', None)
    autoaxis  = kwargs.get('autoaxis', False)
    verbose   = kwargs.get('verbose' , False)
    epssave   = kwargs.get('epssave' , False)
    showplot  = kwargs.get('showplot', False)
    slc_nmb   = kwargs.get('slc_nmb' , None) 
    clp_fnc   = kwargs.get('clp_fnc' ,'sum')

    redshift  = kwargs.get('redshift' ,'1')
    rst_frq   = kwargs.get('rst_frq'  ,'1')

    x_ref    = kwargs.get('x_ref',0)
    y_ref    = kwargs.get('y_ref',0)
    ap_size  = kwargs.get('ap_size',0)


    z_avg     = kwargs.get('z_avg',Header_Get(Cube2bplot_3D,'STZ_AVG'))
    z_med     = kwargs.get('z_med',Header_Get(Cube2bplot_3D,'STZ_MED'))
    frq_r     = kwargs.get('frq_r',Header_Get(Cube2bplot_3D,'RESTFRQ'))
    z_f2l     = z_med

    cube_data     = np.asarray(apgtdt(Cube2bplot_3D,memmap=False) )

    Cube_Info = Cube_Header_Get(Cube2bplot_3D,frq_r* u.Hz)
    FRQ_AXS   = Cube_Info[16].value
    VEL_AXS   = Cube_Info[17].value

    if slc_nmb != None:
        data_2b_plot = cube_data[slc_nmb]
        doublecollapsed = np.ravel(data_2b_plot)
        PlotTitle = 'Slice number: '+ str(slc_nmb+1) + '-' +str(round(VEL_AXS[slc_nmb],0)) + ' km/s'
        PLOTFILENAME = (str(Cube2bplot_3D).split('.fits')[0]).split('/')[-1] + '-3D-slc-'+str(slc_nmb+1)+'.pdf'

        if dest_dir != None:
            PLOTFILENAME = str(dest_dir)  + '/' + PLOTFILENAME
        elif dest_dir == None:
            PLOTFILENAME = plt_dir_res    + '/' + PLOTFILENAME
        Message = 'Generated Plot: ' + PLOTFILENAME + ' slice number : ' + str(slc_nmb+1)

    elif slc_nmb == None:
        PlotTitle = 'Collapse: ' + str(clp_fnc.upper())
        PLOTFILENAME = (str(Cube2bplot_3D).split('.fits')[0]).split('/')[-1] + '-3D-'+str(clp_fnc)+'.pdf'
        if dest_dir != None:
            PLOTFILENAME = str(dest_dir)  + '/' + PLOTFILENAME
        elif dest_dir == None:
            PLOTFILENAME = plt_dir_res    + '/' + PLOTFILENAME
        Message = 'Generated Plot: ' + PLOTFILENAME + ' collapse ('+str(clp_fnc)+')'

        if clp_fnc == 'sum':
            cube_data_clp   = np.asarray(np.nansum(np.array(cube_data)   , axis=0)) 
            cube_data_clp   = cube_data_clp
            data_2b_plot    = cube_data_clp
        elif clp_fnc == 'med':
            cube_data_clp   = np.asarray(np.nanmedian(np.array(cube_data), axis=0)) 
            cube_data_clp   = cube_data_clp
            data_2b_plot    = cube_data_clp
        elif clp_fnc == 'avg':
            cube_data_clp   = np.asarray(np.nanmean(np.array(cube_data)  , axis=0)) 
            cube_data_clp   = cube_data_clp
            data_2b_plot    = cube_data_clp
    freq_num,y_num,x_num = cube_data.shape
    x    = np.arange(0,x_num,1)
    y    = np.arange(0,y_num,1)
    x, y = np.meshgrid(x,y)
    z    = data_2b_plot

    nx_f2DG, ny_f2DG = x_num,y_num 
    nx,ny            = nx_f2DG,ny_f2DG

    X0_f2DG     = kwargs.get('X0_f2DG',int(np.ceil(nx_f2DG/2)))
    Y0_f2DG     = kwargs.get('Y0_f2DG',int(np.ceil(ny_f2DG/2)))

    MAX_CTR     = data_2b_plot[x_num/2,y_num/2]

    fxsize=9
    fysize=8
    f = plt.figure(num=None, figsize=(fxsize, fysize), dpi=180, facecolor='w',
        edgecolor='k')
    plt.subplots_adjust(
        left    = (16/25.4)/fxsize, 
        bottom  = (12/25.4)/fysize, 
        right   = 1 - (6/25.4)/fxsize, 
        top     = 1 - (15/25.4)/fysize)
    plt.subplots_adjust(hspace=0)

    
    ax110 = f.gca(projection='3d')
    ax110.set_rasterization_zorder(1)
    
    plt.title(PlotTitle + ' (' +  str(+x_ref) + ','+str(+y_ref)+') ' + ' $\\varnothing$: '+ str(ap_size) + '"')

    ax110.set_xlabel('X',family = 'serif')
    ax110.set_ylabel('Y',family = 'serif')
    ax110.set_zlabel('S [Jy]',family = 'serif')

    surf = ax110.plot_surface(x,y,z, cmap=cm.viridis,rstride=1, cstride=1,
        linewidth=0, antialiased=False,vmin=np.nanmin(z), vmax=np.nanmax(z))
    divider = make_axes_locatable(ax110)
    cax  = divider.append_axes("right", size="5%", pad=0.05)    
    cbar = plt.colorbar(cax=cax)
    f.colorbar(surf, shrink=0.25, pad=0.05,aspect=5,format=mpl.ticker.FuncFormatter(fmt))

    ax110.xaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-X0_f2DG+x_ref))
    ax110.yaxis.set_major_formatter(ScaledFormatter(dx=1.0,x0=-Y0_f2DG+y_ref))

    ax110.zaxis.set_major_formatter(mpl.ticker.FuncFormatter(fmt))

    plt.savefig(PLOTFILENAME)

    if verbose == True:
        print
        print (colored(Message,'cyan'))
    elif verbose ==False:
        pass
    plt.close('all')

def Plot_Cube_Slices(Cube1,Cube2,Cube3,Cube4,Cube5,Cube6,*args, **kwargs):
    mask_in        = np.asarray((apfts.open(Cube1)[0]).data)
    cub_crc_in     = np.asarray((apfts.open(Cube2)[0]).data)
    cub_crc_wgt_in = np.asarray((apfts.open(Cube3)[0]).data)
    mask_ot        = np.asarray((apfts.open(Cube4)[0]).data)
    cub_crc_ot     = np.asarray((apfts.open(Cube5)[0]).data)
    cub_crc_wgt_ot = np.asarray((apfts.open(Cube6)[0]).data)

    dest_dir_plt   = kwargs.get('dest_dir_plt',None)
    z_avg          = kwargs.get('z_avg',Header_Get(Cube1,'STZ_AVG'))
    z_med          = kwargs.get('z_med',Header_Get(Cube1,'STZ_MED'))
    frq_r          = kwargs.get('frq_r',Header_Get(Cube1,'RESTFRQ'))
    z_f2l          = z_med

    prefix         = kwargs.get('prefix','')

    FLS_in=[Cube1,Cube2,Cube3]
    FLS_in=[Cube4,Cube5,Cube6]

    Cube_Info = Cube_Header_Get(Cube1,frq_r* u.Hz)
    FRQ_AXS   = Cube_Info[16].value
    VEL_AXS   = Cube_Info[17].value

    XAXIS     = VEL_AXS

    if dest_dir_plt != None:
        PLOTFILENAME = dest_dir_plt + prefix + (Cube1.split('.fits',1)[0]).rsplit('/',1)[-1] + '-slices.pdf'
    elif dest_dir_plt == None:
        PLOTFILENAME = ana_dir_plt + prefix + (Cube1.split('.fits',1)[0]).rsplit('/',1)[-1] + '-slices.pdf'

 
    # The PDF document
    pdf_pages         = PdfPages(PLOTFILENAME)
    # Generate the pages
    nb_plots          = mask_in.shape[0]
    nb_plots_per_page = 10
    nb_pages          = int(np.ceil(nb_plots / float(nb_plots_per_page)))

    col_num           = 2
    row_num           = nb_plots_per_page/col_num

    
    X_C = Header_Get(Cube1,'XCT_FIT')
    Y_C = Header_Get(Cube1,'YCT_FIT')
    radii_as_ot= Header_Get(Cube1,'RAD_EXT')
    radii_as_in= Header_Get(Cube5,'RAD_EXT')
    page_num = 0
    for i,(vel,spl1,spl2,spl3,spl4,spl5,spl6) in enumerate(zip(VEL_AXS,mask_in,cub_crc_in,cub_crc_wgt_in,mask_ot,cub_crc_ot,cub_crc_wgt_ot)):
        if i % nb_plots_per_page == 0:
            page_num = page_num +1
            grid_panel=0

            fxsize=22
            fysize=13
            f = plt.figure(num=None, figsize=(fysize, fxsize), dpi=180, facecolor='w',
                edgecolor='k')
            plt.subplots_adjust(
                left    = (22/25.4)/fxsize, 
                bottom  = (6/25.4)/fysize, 
                right   = 1 - (20.4/25.4)/fxsize, 
                top     = 1 - (15/25.4)/fysize)

            if len(cub_crc_in)% 2 == 0:
                col_num = 2
                lin_num = len(cub_crc_in)/2
            else:
                col_num = 2
                lin_num = (len(cub_crc_in)/2) + 1

            outer_grid = gridspec.GridSpec(row_num,col_num,wspace=0.15, hspace=0.15)            
            plt.suptitle((Cube1.split('.fits',1)[0]).rsplit('/',1)[-1] +' '+str(X_C)+','+str(Y_C) +  ' ' +str(radii_as_ot) + 'as ' +str(radii_as_in) + 'as  (' + str(page_num) +')')
        elif (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
            circle=[]
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, 3, 
                subplot_spec=outer_grid[grid_panel],hspace=0.0,wspace=0.0)
            ############################ UNO ##################################
            ax110 = plt.Subplot(f, inner_grid[0,0])

            f.add_subplot(ax110)
            ax110.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax110.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            ax110.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            yticklabels = ax110.get_yticklabels()
            xticklabels = ax110.get_xticklabels()
            plt.setp(yticklabels, visible=True)
            plt.setp(xticklabels, visible=True)
            
            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax110.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')
            plt.ylabel('Y',fontsize=14,family = 'serif')

            circle = plt.imshow(spl4, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube4,'XXT_MIN'), Header_Get(Cube4,'XXT_MAX'),    Header_Get(Cube4,'YXT_MIN'),    Header_Get(Cube4,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='black', marker='x')

            ############################ DOS ##################################
            ax120 = plt.Subplot(f, inner_grid[0,1],sharey=ax110)

            f.add_subplot(ax120)
            ax120.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax120.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='white',bottom=True,top=True,left=True,right=True)
            ax120.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='white',bottom=True,top=True,left=True,right=True)
            ax120.set_title(str(round(VEL_AXS[i],0)) + 'km/s '+str(i+1),fontsize=14,family='serif')
            yticklabels = ax120.get_yticklabels()
            xticklabels = ax120.get_xticklabels()
            plt.setp(yticklabels, visible=False)
            plt.setp(xticklabels, visible=True)     

            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax120.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')

            circle = plt.imshow(spl5, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube5,'XXT_MIN'), Header_Get(Cube5,'XXT_MAX'),    Header_Get(Cube5,'YXT_MIN'),    Header_Get(Cube5,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='white', marker='x')
            ############################ TRES ##################################
            ax130 = plt.Subplot(f, inner_grid[0,2],sharey=ax120)

            f.add_subplot(ax130)
            ax130.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax130.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            ax130.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True,labelleft=False,labelright=True)
            yticklabels = ax130.get_yticklabels()
            xticklabels = ax130.get_xticklabels()
                        
            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax130.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')

            circle = plt.imshow(spl6, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube6,'XXT_MIN'), Header_Get(Cube6,'XXT_MAX'),    Header_Get(Cube6,'YXT_MIN'),    Header_Get(Cube6,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='black', marker='x')
            
            min_y, max_y = ax130.get_ylim()
            min_x, max_x = ax130.get_xlim()

            ############################ UNO ##################################
            ax210 = plt.Subplot(f, inner_grid[1,0],sharex=ax110)

            f.add_subplot(ax210)
            ax210.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax210.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            ax210.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            yticklabels = ax210.get_yticklabels()
            xticklabels = ax210.get_xticklabels()
            
            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax210.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')
            plt.ylabel('Y',fontsize=14,family = 'serif')

            circle = plt.imshow(spl1, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube1,'XXT_MIN'), Header_Get(Cube1,'XXT_MAX'),    Header_Get(Cube1,'YXT_MIN'),    Header_Get(Cube1,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='white', marker='x')

            plt.xlim([min_x,max_x])
            xmin, xmax = plt.xlim()
            plt.xlim((xmin,xmax))

            plt.ylim([min_y,max_y])
            ymin, ymax = plt.ylim()
            plt.ylim((ymin,ymax))

            ############################ DOS ##################################

            ax220 = plt.Subplot(f, inner_grid[1,1],sharey=ax210,sharex=ax120)

            f.add_subplot(ax220)
            ax220.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax220.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            ax220.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            yticklabels = ax220.get_yticklabels()
            xticklabels = ax220.get_xticklabels()
            plt.setp(yticklabels, visible=False)
            plt.setp(xticklabels, visible=True)     

            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax220.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')

            circle = plt.imshow(spl2, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube2,'XXT_MIN'), Header_Get(Cube2,'XXT_MAX'),    Header_Get(Cube2,'YXT_MIN'),    Header_Get(Cube2,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='white', marker='x')

            plt.xlim([min_x,max_x])
            xmin, xmax = plt.xlim()
            plt.xlim((xmin,xmax))

            plt.ylim([min_y,max_y])
            ymin, ymax = plt.ylim()
            plt.ylim((ymin,ymax))

            ############################ TRES ##################################
            ax230 = plt.Subplot(f, inner_grid[1,2],sharey=ax220,sharex=ax130)

            f.add_subplot(ax230)
            ax230.set_rasterization_zorder(1)
            plt.autoscale(enable=True, axis='y', tight=False)
            ax230.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
            ax230.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True,labelleft=False,labelright=True)
            yticklabels = ax230.get_yticklabels()
            xticklabels = ax230.get_xticklabels()
            plt.setp(yticklabels, visible=True)
            plt.setp(xticklabels, visible=True)

            plt.tick_params(which='both', width=1.0)
            plt.tick_params(which='major', length=10,labelsize=5)
            plt.tick_params(which='minor', length=5)
            ax230.minorticks_on()
            plt.xlabel('X',fontsize=14,family = 'serif')

            circle = plt.imshow(spl3, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
                extent = ((Header_Get(Cube3,'XXT_MIN'), Header_Get(Cube3,'XXT_MAX'),    Header_Get(Cube3,'YXT_MIN'),    Header_Get(Cube3,'YXT_MAX'))))
            plt.scatter(X_C, Y_C, s=6, c='white', marker='x')

            plt.xlim([min_x,max_x])
            xmin, xmax = plt.xlim()
            plt.xlim((xmin,xmax))

            plt.ylim([min_y,max_y])
            ymin, ymax = plt.ylim()
            plt.ylim((ymin,ymax))

            ############################ FIN ##################################
            pdf_pages.savefig(f)
            grid_panel=0

        circle=[]
        
        inner_grid = gridspec.GridSpecFromSubplotSpec(2, 3, 
            subplot_spec=outer_grid[grid_panel],hspace=0.0,wspace=0.0)
        
        ############################ UNO ##################################
        ax110 = plt.Subplot(f, inner_grid[0,0])

        f.add_subplot(ax110)
        ax110.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax110.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax110.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        yticklabels = ax110.get_yticklabels()
        xticklabels = ax110.get_xticklabels()
        plt.setp(yticklabels, visible=True)
        plt.setp(xticklabels, visible=True)
        
        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax110.minorticks_on()
        plt.xlabel('X',fontsize=14,family = 'serif')
        plt.ylabel('Y',fontsize=14,family = 'serif')

        circle = plt.imshow(spl4, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube4,'XXT_MIN'), Header_Get(Cube4,'XXT_MAX'),    Header_Get(Cube4,'YXT_MIN'),    Header_Get(Cube4,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='black', marker='x')
        ############################ DOS ##################################

        ax120 = plt.Subplot(f, inner_grid[0,1],sharey=ax110)

        f.add_subplot(ax120)
        ax120.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax120.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='white',bottom=True,top=True,left=True,right=True)
        ax120.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='white',bottom=True,top=True,left=True,right=True)
        ax120.set_title(str(round(VEL_AXS[i],0)) + 'km/s '+str(i+1),fontsize=14,family='serif')
        yticklabels = ax120.get_yticklabels()
        xticklabels = ax120.get_xticklabels()
        plt.setp(yticklabels, visible=False)
        plt.setp(xticklabels, visible=True)     

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax120.minorticks_on()
        plt.xlabel('X',fontsize=14,family = 'serif')

        circle = plt.imshow(spl5, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube5,'XXT_MIN'), Header_Get(Cube5,'XXT_MAX'),    Header_Get(Cube5,'YXT_MIN'),    Header_Get(Cube5,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='white', marker='x')
        ############################ TRES ##################################
        ax130 = plt.Subplot(f, inner_grid[0,2],sharey=ax120)

        f.add_subplot(ax130)
        ax130.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax130.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax130.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True,labelleft=False,labelright=True)
        yticklabels = ax130.get_yticklabels()
        xticklabels = ax130.get_xticklabels()
                    
        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax130.minorticks_on()
        plt.xlabel('X',fontsize=14,family = 'serif')

        circle = plt.imshow(spl6, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube6,'XXT_MIN'), Header_Get(Cube6,'XXT_MAX'),    Header_Get(Cube6,'YXT_MIN'),    Header_Get(Cube6,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='black', marker='x')
        
        min_y, max_y = ax130.get_ylim()
        min_x, max_x = ax130.get_xlim()

        ############################ UNO ##################################
        ax210 = plt.Subplot(f, inner_grid[1,0],sharex=ax110)

        f.add_subplot(ax210)
        ax210.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax210.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax210.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        yticklabels = ax210.get_yticklabels()
        xticklabels = ax210.get_xticklabels()
        
        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax210.minorticks_on()
        plt.xlabel('X',fontsize=14,family = 'serif')
        plt.ylabel('Y',fontsize=14,family = 'serif')

        circle = plt.imshow(spl1, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube1,'XXT_MIN'), Header_Get(Cube1,'XXT_MAX'),    Header_Get(Cube1,'YXT_MIN'),    Header_Get(Cube1,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='white', marker='x')

        plt.xlim([min_x,max_x])
        xmin, xmax = plt.xlim()
        plt.xlim((xmin,xmax))

        plt.ylim([min_y,max_y])
        ymin, ymax = plt.ylim()
        plt.ylim((ymin,ymax))

        ############################ DOS ##################################

        ax220 = plt.Subplot(f, inner_grid[1,1],sharey=ax210,sharex=ax120)

        f.add_subplot(ax220)
        ax220.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax220.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax220.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        #ax220.set_title(str(slice_number+1),family='serif')
        yticklabels = ax220.get_yticklabels()
        xticklabels = ax220.get_xticklabels()
        plt.setp(yticklabels, visible=False)
        plt.setp(xticklabels, visible=True)     

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax220.minorticks_on()

        circle = plt.imshow(spl2, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube2,'XXT_MIN'), Header_Get(Cube2,'XXT_MAX'),    Header_Get(Cube2,'YXT_MIN'),    Header_Get(Cube2,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='white', marker='x')
        
        plt.xlim([min_x,max_x])
        xmin, xmax = plt.xlim()
        plt.xlim((xmin,xmax))

        plt.ylim([min_y,max_y])
        ymin, ymax = plt.ylim()
        plt.ylim((ymin,ymax))

        ############################ TRES ##################################
        ax230 = plt.Subplot(f, inner_grid[1,2],sharey=ax220,sharex=ax130)

        f.add_subplot(ax230)
        ax230.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax230.xaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax230.yaxis.set_tick_params(which='both',labelsize=4,direction='in',color='black',bottom=True,top=True,left=True,right=True,labelleft=False,labelright=True)
        yticklabels = ax230.get_yticklabels()
        xticklabels = ax230.get_xticklabels()
        plt.setp(yticklabels, visible=True)
        plt.setp(xticklabels, visible=True)

                    
        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10,labelsize=5)
        plt.tick_params(which='minor', length=5)
        ax230.minorticks_on()
        plt.xlabel('X',fontsize=14,family = 'serif')

        circle = plt.imshow(spl3, cmap=plt.cm.viridis,interpolation='nearest', origin='lower',
            extent = ((Header_Get(Cube3,'XXT_MIN'), Header_Get(Cube3,'XXT_MAX'),    Header_Get(Cube3,'YXT_MIN'),    Header_Get(Cube3,'YXT_MAX'))))
        plt.scatter(X_C, Y_C, s=6, c='white', marker='x')


        plt.xlim([min_x,max_x])
        xmin, xmax = plt.xlim()
        plt.xlim((xmin,xmax))

        plt.ylim([min_y,max_y])
        ymin, ymax = plt.ylim()
        plt.ylim((ymin,ymax))

        ############################ FIN ##################################
        grid_panel=grid_panel+1

    pdf_pages.close()
    print (colored('Circular slices on: '+PLOTFILENAME,'cyan'))

def Cube_Header_Get(cube_header_ipt,freq_rfr,*args, **kwargs):
    verbose    = kwargs.get('verbose',False) 
    redshift   = kwargs.get('redshift',0)
    freq_step  = kwargs.get('freq_step',Header_Get(cube_header_ipt,'CDELT3'))
    freq_init  = kwargs.get('freq_init',Header_Get(cube_header_ipt,'CRVAL3'))
    freq_obs_f = kwargs.get('freq_obs_f',Header_Get(cube_header_ipt,'RESTFRQ'))#freq_obs
    freq_obs   = kwargs.get('freq_obs',Redshifted_freq(freq_rfr,redshift))    #freq_obs_f 

    vel_step   = kwargs.get('vel_step',None)

    DIM1 = Header_Get(cube_header_ipt,'NAXIS1')
    DIM2 = Header_Get(cube_header_ipt,'NAXIS2')
    DIM3 = Header_Get(cube_header_ipt,'NAXIS3')

    RVL1 = Header_Get(cube_header_ipt,'CRVAL1')
    RVL2 = Header_Get(cube_header_ipt,'CRVAL2')
    RVL3 = Header_Get(cube_header_ipt,'CRVAL3')

    RES1 = Header_Get(cube_header_ipt,'CDELT1')
    RES2 = Header_Get(cube_header_ipt,'CDELT2')

    cube                      = scspc.read(cube_header_ipt) 
    freq_step_shifted         = Redshifted_freq(freq_step,float(redshift))
    freq_init_shifted         = Redshifted_freq(freq_init,float(redshift))
    
    vel_step                  = abs(Thermal_Dopp_vel(freq_obs_f,freq_obs_f-abs((freq_step)),redshift_freq=redshift))
    vel_step_frq_bck          = Thermal_Dopp_freq(freq_obs_f,vel_step,redshift_freq=redshift) - freq_obs_f

    cube                      = cube.with_spectral_unit(u.Hz,velocity_convention='radio',rest_value=freq_obs_f* u.Hz)
    cube_vel                  = cube.with_spectral_unit(u.km / u.s,velocity_convention='radio')

    freq_max,freq_min = cube.spectral_extrema[-1].value,cube.spectral_extrema[0].value
    vel_max,vel_min   = cube_vel.spectral_extrema[-1].value,cube_vel.spectral_extrema[0].value
    
    if verbose == True:
        print (colored('Reading cube: ' + cube_header_ipt,'yellow'))
        print
        print ('Cube                                         : ',cube_header_ipt)
        print
        print (colored(cube,'yellow'))
        print
        print ('Dimensions                                   : ',DIM1,'X',DIM2,'X',DIM3)
        print ('CII restframe                                : ',freq_rfr)
        print ('CII will be observed at (fits file)          : ',freq_obs_f)
        print ('CII will be observed at                      : ',freq_obs)
        print
        print ('Frequency step                               : ',freq_step)
        print ('Initial frequency                            : ',freq_init)
        print ('Max,Min frequencies                          : ',freq_max,freq_min)
        print ('Velocity  step from Frequency        [kms-1] : ',vel_step)
        print ('Frequency step from Velocity                 : ',vel_step_frq_bck)
        print ('Max,Min velocities                   [kms-1] : ',vel_max,vel_min)
    elif verbose == False:
        pass
    return freq_rfr,freq_obs_f,DIM1,DIM2,DIM3,RVL1,RVL2,RVL3,RES1,RES2,freq_step,vel_step,freq_max,freq_min,vel_max,vel_min,cube.spectral_axis,cube_vel.spectral_axis

def MCMC_Confidence_Interval(iterations_mc,line1,line2,method,error,z_1,v_1,z_2,v_2,var_smpls_mcmc,nmb_smpls_mcmc,z_lne1_mcmc,z_lne2_mcmc,v_lne1_mcmc,v_lne2_mcmc,spec_file_plt_lne1_avg_mcmc,spec_file_plt_lne1_med_mcmc,spec_file_plt_lne2_avg_mcmc,spec_file_plt_lne2_med_mcmc,flx_lne1_avg_hdr,flx_lne1_med_hdr,flx_lne1_avg_hdr_e,flx_lne1_med_hdr_e,flx_lne2_avg_hdr,flx_lne2_med_hdr,flx_lne2_avg_hdr_e,flx_lne2_med_hdr_e,*args, **kwargs):
    dest_dir_plt = kwargs.get('dest_dir_plt',None)
    nbins        = kwargs.get('nbins' ,200)
    sfx_lne1     = kwargs.get('sfx_lne1' ,'12CO')
    sfx_lne2     = kwargs.get('sfx_lne2' ,'13CO')
    label1       = kwargs.get('label1','12CO')
    label2       = kwargs.get('label2','13CO')
    label3       = kwargs.get('label2','LFIR')
    line1        = kwargs.get('line1','12CO')
    line2        = kwargs.get('line2','13CO')

    plt_scl      = kwargs.get('plt_scl',None)
    log_lm       = kwargs.get('log_lm','error')
    func1        = kwargs.get('func1','avg')
    func2        = kwargs.get('func2','med')

    flx_hdr     = kwargs.get('flx_hdr',None)
    lum_hdr     = kwargs.get('lum_hdr',None)
    flx_hdr_e   = kwargs.get('flx_hdr_e',None)
    lum_hdr_e1  = kwargs.get('lum_hdr_e1',None)
    lum_hdr_e2  = kwargs.get('lum_hdr_e2',None)
    lum_hdr_e11 = kwargs.get('lum_hdr_e11',None)
    lum_hdr_e21 = kwargs.get('lum_hdr_e21',None)
    lum_hdr_e12 = kwargs.get('lum_hdr_e12',None)
    lum_hdr_e22 = kwargs.get('lum_hdr_e22',None)
    lum_hdr_e13 = kwargs.get('lum_hdr_e13',None)
    lum_hdr_e23 = kwargs.get('lum_hdr_e23',None)

    plot_dist_hist = kwargs.get('plot_dist_hist',True)
    if line1 == '13CO':
        restframe_frequency_1      =   110.20137E9           
    elif line1 == '12CO':
        restframe_frequency_1      =   115.271208E9
    elif line1 == '18CO':
        restframe_frequency_1      =   109.78217340E9

    if line2 == '13CO':
        restframe_frequency_2      =   110.20137E9           
    elif line2 == '12CO':
        restframe_frequency_2      =   115.271208E9
    elif line2 == '18CO':
        restframe_frequency_2      =   109.78217340E9

    Splt_Vars_Plt = split_variable_vars(var_smpls_mcmc)
    Splt_Col = Splt_Vars_Plt[0]
    Splt_Hdr = Splt_Vars_Plt[1]
    Splt_Hdr_Cmt = Splt_Vars_Plt[2]
    Splt_CNm = Splt_Vars_Plt[3]
    Splt_Hdr_Plt = Splt_Vars_Plt[4]
    Splt_Plt_lbl = Splt_Vars_Plt[5]

    rds_hdr     = 'STZ_MED'
    vrx_hdr     = 'STS_MED'

    ############ARRAYS FOR TABLES FOR MCMC################
    Z1_AVG_E1_MC   , Z1_AVG_E2_MC   = [], []
    Z2_AVG_E1_MC   , Z2_AVG_E2_MC   = [], []
    Z1_MED_E1_MC   , Z1_MED_E2_MC   = [], []
    Z2_MED_E1_MC   , Z2_MED_E2_MC   = [], []
    V1_AVG_E1_MC   , V1_AVG_E2_MC   = [], []
    V2_AVG_E1_MC   , V2_AVG_E2_MC   = [], []
    V1_MED_E1_MC   , V1_MED_E2_MC   = [], []
    V2_MED_E1_MC   , V2_MED_E2_MC   = [], []

    S1_AVG_E1_MC   ,  S2_AVG_E1_MC   = [], []
    S1_AVG_E2_MC   ,  S2_AVG_E2_MC   = [], []
    S1_MED_E1_MC   ,  S2_MED_E1_MC   = [], []
    S1_MED_E2_MC   ,  S2_MED_E2_MC   = [], []

    L1_AVG_E1_MC_1 , L2_AVG_E1_MC_1  = [], []
    L1_MED_E1_MC_1 , L2_MED_E1_MC_1  = [], []
    L1_AVG_E1_MC_2 , L2_AVG_E1_MC_2  = [], []
    L1_MED_E1_MC_2 , L2_MED_E1_MC_2  = [], []

    L1_AVG_E2_MC_1 , L2_AVG_E2_MC_1  = [], []
    L1_MED_E2_MC_1 , L2_MED_E2_MC_1  = [], []
    L1_AVG_E2_MC_2 , L2_AVG_E2_MC_2  = [], []
    L1_MED_E2_MC_2 , L2_MED_E2_MC_2  = [], []

    S12_AVG_E1_MC  , S12_AVG_E2_MC   = [], []
    S12_MED_E1_MC  , S12_MED_E2_MC   = [], []

    L12_AVG_E1_MC_1, L12_AVG_E2_MC_1 = [], []
    L12_MED_E1_MC_1, L12_MED_E2_MC_1 = [], []
    L12_AVG_E1_MC_2, L12_AVG_E2_MC_2 = [], []
    L12_MED_E1_MC_2, L12_MED_E2_MC_2 = [], []
    ###################TABLES########################
    SMPL = []

    FA1_MSR                                  = []
    FA1_MSR_MC_1_SG_L, FA1_MSR_MC_1_SG_H     = [], []
    FA1_MSR_MC_2_SG_L, FA1_MSR_MC_2_SG_H     = [], []
    FA1_MSR_MC_3_SG_L, FA1_MSR_MC_3_SG_H     = [], []

    FA2_MSR                                  = []
    FA2_MSR_MC_1_SG_L, FA2_MSR_MC_1_SG_H     = [], []
    FA2_MSR_MC_2_SG_L, FA2_MSR_MC_2_SG_H     = [], []
    FA2_MSR_MC_3_SG_L, FA2_MSR_MC_3_SG_H     = [], []

    FA12_MSR                                 = []
    FA12_MSR_MC_1_SG_L, FA12_MSR_MC_1_SG_H   = [], []
    FA12_MSR_MC_2_SG_L, FA12_MSR_MC_2_SG_H   = [], []
    FA12_MSR_MC_3_SG_L, FA12_MSR_MC_3_SG_H   = [], []

    FM1_MSR                                  = []
    FM1_MSR_MC_1_SG_L, FM1_MSR_MC_1_SG_H     = [], []
    FM1_MSR_MC_2_SG_L, FM1_MSR_MC_2_SG_H     = [], []
    FM1_MSR_MC_3_SG_L, FM1_MSR_MC_3_SG_H     = [], []

    FM2_MSR                                  = []
    FM2_MSR_MC_1_SG_L, FM2_MSR_MC_1_SG_H     = [], []
    FM2_MSR_MC_2_SG_L, FM2_MSR_MC_2_SG_H     = [], []
    FM2_MSR_MC_3_SG_L, FM2_MSR_MC_3_SG_H     = [], []

    FM12_MSR                                 = []
    FM12_MSR_MC_1_SG_L, FM12_MSR_MC_1_SG_H   = [], []
    FM12_MSR_MC_2_SG_L, FM12_MSR_MC_2_SG_H   = [], []
    FM12_MSR_MC_3_SG_L, FM12_MSR_MC_3_SG_H   = [], []

    LA1_MSR                                  = []
    LA1_MSR_MC_1_SG_L, LA1_MSR_MC_1_SG_H     = [], []
    LA1_MSR_MC_2_SG_L, LA1_MSR_MC_2_SG_H     = [], []
    LA1_MSR_MC_3_SG_L, LA1_MSR_MC_3_SG_H     = [], []

    LA2_MSR                                  = []
    LA2_MSR_MC_1_SG_L, LA2_MSR_MC_1_SG_H     = [], []
    LA2_MSR_MC_2_SG_L, LA2_MSR_MC_2_SG_H     = [], []
    LA2_MSR_MC_3_SG_L, LA2_MSR_MC_3_SG_H     = [], []

    LA12_MSR                                 = []
    LA12_MSR_MC_1_SG_L, LA12_MSR_MC_1_SG_H   = [], []
    LA12_MSR_MC_2_SG_L, LA12_MSR_MC_2_SG_H   = [], []
    LA12_MSR_MC_3_SG_L, LA12_MSR_MC_3_SG_H   = [], []

    LM1_MSR                                  = []
    LM1_MSR_MC_1_SG_L, LM1_MSR_MC_1_SG_H     = [], []
    LM1_MSR_MC_2_SG_L, LM1_MSR_MC_2_SG_H     = [], []
    LM1_MSR_MC_3_SG_L, LM1_MSR_MC_3_SG_H     = [], []

    LM2_MSR                                  = []
    LM2_MSR_MC_1_SG_L, LM2_MSR_MC_1_SG_H     = [], []
    LM2_MSR_MC_2_SG_L, LM2_MSR_MC_2_SG_H     = [], []
    LM2_MSR_MC_3_SG_L, LM2_MSR_MC_3_SG_H     = [], []

    LM12_MSR                                 = []
    LM12_MSR_MC_1_SG_L, LM12_MSR_MC_1_SG_H   = [], []
    LM12_MSR_MC_2_SG_L, LM12_MSR_MC_2_SG_H   = [], []
    LM12_MSR_MC_3_SG_L, LM12_MSR_MC_3_SG_H   = [], []

    LLA1_MSR                                 = []
    LLA1_MSR_MC_1_SG_L, LLA1_MSR_MC_1_SG_H   = [], []
    LLA1_MSR_MC_2_SG_L, LLA1_MSR_MC_2_SG_H   = [], []
    LLA1_MSR_MC_3_SG_L, LLA1_MSR_MC_3_SG_H   = [], []

    LLA2_MSR                                 = []
    LLA2_MSR_MC_1_SG_L, LLA2_MSR_MC_1_SG_H   = [], []
    LLA2_MSR_MC_2_SG_L, LLA2_MSR_MC_2_SG_H   = [], []
    LLA2_MSR_MC_3_SG_L, LLA2_MSR_MC_3_SG_H   = [], []

    LLA12_MSR                                = []
    LLA12_MSR_MC_1_SG_L, LLA12_MSR_MC_1_SG_H = [], []
    LLA12_MSR_MC_2_SG_L, LLA12_MSR_MC_2_SG_H = [], []
    LLA12_MSR_MC_3_SG_L, LLA12_MSR_MC_3_SG_H = [], []

    LLM1_MSR                                 = []
    LLM1_MSR_MC_1_SG_L, LLM1_MSR_MC_1_SG_H   = [], []
    LLM1_MSR_MC_2_SG_L, LLM1_MSR_MC_2_SG_H   = [], []
    LLM1_MSR_MC_3_SG_L, LLM1_MSR_MC_3_SG_H   = [], []

    LLM2_MSR                                 = []
    LLM2_MSR_MC_1_SG_L, LLM2_MSR_MC_1_SG_H   = [], []
    LLM2_MSR_MC_2_SG_L, LLM2_MSR_MC_2_SG_H   = [], []
    LLM2_MSR_MC_3_SG_L, LLM2_MSR_MC_3_SG_H   = [], []

    LLM12_MSR                                = []
    LLM12_MSR_MC_1_SG_L, LLM12_MSR_MC_1_SG_H = [], []
    LLM12_MSR_MC_2_SG_L, LLM12_MSR_MC_2_SG_H = [], []
    LLM12_MSR_MC_3_SG_L, LLM12_MSR_MC_3_SG_H = [], []

    ZA1_MSR                                  = []
    ZA1_MSR_MC_1_SG_L, ZA1_MSR_MC_1_SG_H     = [], []
    ZA1_MSR_MC_2_SG_L, ZA1_MSR_MC_2_SG_H     = [], []
    ZA1_MSR_MC_3_SG_L, ZA1_MSR_MC_3_SG_H     = [], []

    ZA2_MSR                                  = []
    ZA2_MSR_MC_1_SG_L, ZA2_MSR_MC_1_SG_H     = [], []
    ZA2_MSR_MC_2_SG_L, ZA2_MSR_MC_2_SG_H     = [], []
    ZA2_MSR_MC_3_SG_L, ZA2_MSR_MC_3_SG_H     = [], []

    ZA12_MSR                                 = []
    ZA12_MSR_MC_1_SG_L, ZA12_MSR_MC_1_SG_H   = [], []
    ZA12_MSR_MC_2_SG_L, ZA12_MSR_MC_2_SG_H   = [], []
    ZA12_MSR_MC_3_SG_L, ZA12_MSR_MC_3_SG_H   = [], []

    ZM1_MSR                                  = []
    ZM1_MSR_MC_1_SG_L, ZM1_MSR_MC_1_SG_H     = [], []
    ZM1_MSR_MC_2_SG_L, ZM1_MSR_MC_2_SG_H     = [], []
    ZM1_MSR_MC_3_SG_L, ZM1_MSR_MC_3_SG_H     = [], []

    ZM2_MSR                                  = []
    ZM2_MSR_MC_1_SG_L, ZM2_MSR_MC_1_SG_H     = [], []
    ZM2_MSR_MC_2_SG_L, ZM2_MSR_MC_2_SG_H     = [], []
    ZM2_MSR_MC_3_SG_L, ZM2_MSR_MC_3_SG_H     = [], []

    ZM12_MSR                                 = []
    ZM12_MSR_MC_1_SG_L, ZM12_MSR_MC_1_SG_H   = [], []
    ZM12_MSR_MC_2_SG_L, ZM12_MSR_MC_2_SG_H   = [], []
    ZM12_MSR_MC_3_SG_L, ZM12_MSR_MC_3_SG_H   = [], []

    VA1_MSR                                  = []
    VA1_MSR_MC_1_SG_L, VA1_MSR_MC_1_SG_H     = [], []
    VA1_MSR_MC_2_SG_L, VA1_MSR_MC_2_SG_H     = [], []
    VA1_MSR_MC_3_SG_L, VA1_MSR_MC_3_SG_H     = [], []

    VA2_MSR                                  = []
    VA2_MSR_MC_1_SG_L, VA2_MSR_MC_1_SG_H     = [], []
    VA2_MSR_MC_2_SG_L, VA2_MSR_MC_2_SG_H     = [], []
    VA2_MSR_MC_3_SG_L, VA2_MSR_MC_3_SG_H     = [], []

    VA12_MSR                                 = []
    VA12_MSR_MC_1_SG_L, VA12_MSR_MC_1_SG_H   = [], []
    VA12_MSR_MC_2_SG_L, VA12_MSR_MC_2_SG_H   = [], []
    VA12_MSR_MC_3_SG_L, VA12_MSR_MC_3_SG_H   = [], []

    VM1_MSR                                  = []
    VM1_MSR_MC_1_SG_L, VM1_MSR_MC_1_SG_H     = [], []
    VM1_MSR_MC_2_SG_L, VM1_MSR_MC_2_SG_H     = [], []
    VM1_MSR_MC_3_SG_L, VM1_MSR_MC_3_SG_H     = [], []

    VM2_MSR                                  = []
    VM2_MSR_MC_1_SG_L, VM2_MSR_MC_1_SG_H     = [], []
    VM2_MSR_MC_2_SG_L, VM2_MSR_MC_2_SG_H     = [], []
    VM2_MSR_MC_3_SG_L, VM2_MSR_MC_3_SG_H     = [], []

    VM12_MSR                                 = []
    VM12_MSR_MC_1_SG_L, VM12_MSR_MC_1_SG_H   = [], []
    VM12_MSR_MC_2_SG_L, VM12_MSR_MC_2_SG_H   = [], []
    VM12_MSR_MC_3_SG_L, VM12_MSR_MC_3_SG_H   = [], []
    ###################TABLES########################
    ############ARRAYS FOR TABLES FOR MCMC################

    ############################################ MC ############################################
    print (colored('Confidence intervals through MC','yellow'))
    import random
    n_bins = nbins 

    ###########################################MC DISTRIBUTION ###########################################
    n_z_lne1, bins_z_lne1, patches_z_lne1 = plt.hist(z_lne1_mcmc, n_bins, density=1, histtype='step',cumulative=True)
    n_z_lne2, bins_z_lne2, patches_z_lne2 = plt.hist(z_lne2_mcmc, n_bins, density=1, histtype='step',cumulative=True)
    n_v_lne1, bins_v_lne1, patches_v_lne1 = plt.hist(v_lne1_mcmc, n_bins, density=1, histtype='step',cumulative=True)
    n_v_lne2, bins_v_lne2, patches_v_lne2 = plt.hist(v_lne2_mcmc, n_bins, density=1, histtype='step',cumulative=True)

    p_z_lne1,p_z_lne2,p_v_lne1,p_v_lne2=[],[],[],[]
    for j in range(iterations_mc):
        p_z_lne1.append(random.random())
        p_z_lne2.append(random.random())
        p_v_lne1.append(random.random())
        p_v_lne2.append(random.random())

    index_z_lne1 = np.digitize(p_z_lne1,n_z_lne1)
    index_z_lne2 = np.digitize(p_z_lne2,n_z_lne2)
    index_v_lne2 = np.digitize(p_v_lne1,n_v_lne1)
    index_v_lne2 = np.digitize(p_v_lne2,n_v_lne2)

    i_z_lne1_mc  = index_z_lne1          
    p_z_lne1_mc  = n_z_lne1[index_z_lne1]     
    z_z_lne1_mc  = bins_z_lne1[index_z_lne1]  

    i_z_lne2_mc  = index_z_lne2          
    p_z_lne2_mc  = n_z_lne2[index_z_lne2]     
    z_z_lne2_mc  = bins_z_lne2[index_z_lne2]  


    i_v_lne2_mc  = index_v_lne2          
    p_v_lne1_mc  = n_z_lne1[index_v_lne2]     
    v_z_lne1_mc  = bins_v_lne2[index_v_lne2]  

    i_v_lne2_mc  = index_v_lne2          
    p_v_lne2_mc  = n_z_lne2[index_v_lne2]     
    v_z_lne2_mc  = bins_v_lne2[index_v_lne2]  

    z_lne12_avg_mc = z_z_lne1_mc / z_z_lne2_mc
    z_lne12_med_mc = z_z_lne1_mc / z_z_lne2_mc

    v_lne12_avg_mc = v_z_lne1_mc / v_z_lne2_mc
    v_lne12_med_mc = v_z_lne1_mc / v_z_lne2_mc

    mu_lne1_avg  = Header_Get(spec_file_plt_lne1_avg_mcmc,flx_lne1_avg_hdr)
    mu_lne1_med  = Header_Get(spec_file_plt_lne1_med_mcmc,flx_lne1_med_hdr)
    sgm_lne1_avg = Header_Get(spec_file_plt_lne1_avg_mcmc,flx_lne1_avg_hdr_e)
    sgm_lne1_med = Header_Get(spec_file_plt_lne1_med_mcmc,flx_lne1_med_hdr_e)

    mu_lne2_avg  = Header_Get(spec_file_plt_lne2_avg_mcmc,flx_lne2_avg_hdr)
    mu_lne2_med  = Header_Get(spec_file_plt_lne2_med_mcmc,flx_lne2_med_hdr)
    sgm_lne2_avg = Header_Get(spec_file_plt_lne2_avg_mcmc,flx_lne2_avg_hdr_e)
    sgm_lne2_med = Header_Get(spec_file_plt_lne2_med_mcmc,flx_lne2_med_hdr_e)

    #########################################################
    ##########1DGF WITH NANs DELTA LIKE STD STD_ERR##########
    if np.isnan(sgm_lne1_avg)==True:
        print
        print ('sgm_lne1_avg: ',sgm_lne1_avg)
        print (colored('NaN encountered in sgm value from fits header: ' + str(flx_lne1_avg_hdr_e) + ' for generating de gaussian distribution!','yellow'))
        print (colored(spec_file_plt_lne1_avg_mcmc,'cyan'))
        print (colored('Using alternative error estimation: FTS_AAE.','yellow'))
        print
        sgm_lne1_avg = element[0]
        sgm_lne1_avg = Header_Get(spec_file_plt_lne1_avg_mcmc,'FTS_AAE')

    elif np.isnan(sgm_lne1_avg)==False:
        pass
    if np.isnan(sgm_lne1_med)==True:
        print
        print ('sgm_lne1_med: ',sgm_lne1_med)
        print (colored('NaN encountered in sgm value from fits header: ' + str(flx_lne1_med_hdr_e) + ' for generating de gaussian distribution!','yellow'))
        print (colored(spec_file_plt_lne1_med_mcmc,'cyan'))
        print (colored('Using alternative error estimation: FTS_AAE.','yellow'))
        print
        sgm_lne1_med = element[0]
        sgm_lne1_med = Header_Get(spec_file_plt_lne1_med_mcmc,'FTS_AAE')
    elif np.isnan(sgm_lne1_med)==False:
        pass
    if np.isnan(sgm_lne2_avg)==True:
        print
        print ('sgm_lne2_avg: ',sgm_lne2_avg)
        print (colored('NaN encountered in sgm value from fits header: ' + str(flx_lne2_avg_hdr_e) + ' for generating de gaussian distribution!','yellow'))
        print (colored(spec_file_plt_lne2_avg_mcmc,'cyan'))
        print (colored('Using alternative error estimation: FTS_AAE.','yellow'))
        print
        sgm_lne2_avg = element[0]
        sgm_lne2_avg = Header_Get(spec_file_plt_lne2_avg_mcmc,'FTS_AAE')
    elif np.isnan(sgm_lne2_avg)==False:
        pass
    if np.isnan(sgm_lne2_med)==True:
        print
        print ('sgm_lne2_med: ',sgm_lne2_med)
        print (colored('NaN encountered in sgm value from fits header: ' + str(flx_lne2_med_hdr_e) + ' for generating de gaussian distribution!','yellow'))
        print (colored(spec_file_plt_lne2_med_mcmc,'cyan'))
        print (colored('Using alternative error estimation: FTS_AAE.','yellow'))
        print
        sgm_lne2_med = element[0]
        sgm_lne2_med = Header_Get(spec_file_plt_lne2_med_mcmc,'FTS_AAE')
    elif np.isnan(sgm_lne2_med)==False:
        pass            
    #########################################################
    s_lne1_avg_mc = np.random.normal(mu_lne1_avg, abs(sgm_lne1_avg), iterations_mc)
    s_lne1_med_mc = np.random.normal(mu_lne1_med, abs(sgm_lne1_med), iterations_mc)
    s_lne2_avg_mc = np.random.normal(mu_lne2_avg, abs(sgm_lne2_avg), iterations_mc)
    s_lne2_med_mc = np.random.normal(mu_lne2_med, abs(sgm_lne2_med), iterations_mc)

    try:
        s_lne12_avg_mc = s_lne1_avg_mc / s_lne2_avg_mc
    except RuntimeWarning:
        print ('Error')
        print (s_lne1_avg_mc / s_lne2_avg_mc)
        print
        print (spec_file_plt_lne1_avg_mcmc)
        print (spec_file_plt_lne2_avg_mcmc)
        print (flx_lne1_avg_hdr,flx_lne2_avg_hdr)
        print (flx_lne1_avg_hdr_e,flx_lne2_avg_hdr_e)
        print (mu_lne1_avg)
        print (sgm_lne1_avg)
        print (mu_lne2_avg)
        print (sgm_lne2_avg)
        print
        print ('Quitting!')
        quit()
    try:
        s_lne12_med_mc = s_lne1_med_mc / s_lne2_med_mc
    except RuntimeWarning:
        print ('Error')
        print (s_lne1_med_mc / s_lne2_med_mc)
        print
        print (spec_file_plt_lne1_avg_mcmc)
        print (spec_file_plt_lne2_avg_mcmc)
        print (flx_lne1_med_hdr,flx_lne2_med_hdr)
        print (flx_lne1_med_hdr_e,flx_lne2_med_hdr_e)
        print (mu_lne1_med)
        print (sgm_lne1_med)
        print (mu_lne2_med)
        print (sgm_lne2_med)
        print
        print ('Quitting!'       )
        quit()

    Lum_lne1_avg_mc,Lum_lne1_med_mc = [], []
    Lum_lne2_avg_mc,Lum_lne2_med_mc = [], []
    [Lum_lne1_avg_mc.append(FluxToLum(s_rep,z_rep,restframe_frequency_1)) for s_rep,z_rep in zip(s_lne1_avg_mc,z_z_lne1_mc)]
    [Lum_lne1_med_mc.append(FluxToLum(s_rep,z_rep,restframe_frequency_1)) for s_rep,z_rep in zip(s_lne1_med_mc,z_z_lne1_mc)]
    [Lum_lne2_avg_mc.append(FluxToLum(s_rep,z_rep,restframe_frequency_2)) for s_rep,z_rep in zip(s_lne2_avg_mc,z_z_lne2_mc)]
    [Lum_lne2_med_mc.append(FluxToLum(s_rep,z_rep,restframe_frequency_2)) for s_rep,z_rep in zip(s_lne2_med_mc,z_z_lne2_mc)]

    Lum_lne1_avg_mc = np.asarray(Lum_lne1_avg_mc)
    Lum_lne1_med_mc = np.asarray(Lum_lne1_med_mc)
    Lum_lne2_avg_mc = np.asarray(Lum_lne2_avg_mc)
    Lum_lne2_med_mc = np.asarray(Lum_lne2_med_mc)

    Lum_lne1_avg_mc_1 = Lum_lne1_avg_mc[:,0]
    Lum_lne1_avg_mc_2 = Lum_lne1_avg_mc[:,1]
    Lum_lne1_med_mc_1 = Lum_lne1_med_mc[:,0]
    Lum_lne1_med_mc_2 = Lum_lne1_med_mc[:,1]

    Lum_lne2_avg_mc_1 = Lum_lne2_avg_mc[:,0]
    Lum_lne2_avg_mc_2 = Lum_lne2_avg_mc[:,1]
    Lum_lne2_med_mc_1 = Lum_lne2_med_mc[:,0]
    Lum_lne2_med_mc_2 = Lum_lne2_med_mc[:,1]

    if error ==1:
        e_l_mc = 15.9
        e_s_mc = 84.1
    elif error ==2:
        e_l_mc = 2.30
        e_s_mc = 97.7
    elif error == 3:    
        e_l_mc = 0.20
        e_s_mc = 99.8

    Lum_lne1_avg_mc_1_indx = np.where(Lum_lne1_avg_mc_1==88888.0)
    Lum_lne1_med_mc_1_indx = np.where(Lum_lne1_med_mc_1==88888.0)
    Lum_lne1_avg_mc_2_indx = np.where(Lum_lne1_avg_mc_2==88888.0)
    Lum_lne1_med_mc_2_indx = np.where(Lum_lne1_med_mc_2==88888.0)
    Lum_lne2_avg_mc_1_indx = np.where(Lum_lne2_avg_mc_1==88888.0)
    Lum_lne2_med_mc_1_indx = np.where(Lum_lne2_med_mc_1==88888.0)
    Lum_lne2_avg_mc_2_indx = np.where(Lum_lne2_avg_mc_2==88888.0)
    Lum_lne2_med_mc_2_indx = np.where(Lum_lne2_med_mc_2==88888.0)

    Lum_lne1_avg_mc_1[Lum_lne1_avg_mc_1_indx] = np.nan
    Lum_lne1_med_mc_1[Lum_lne1_med_mc_1_indx] = np.nan
    Lum_lne1_avg_mc_2[Lum_lne1_avg_mc_2_indx] = np.nan
    Lum_lne1_med_mc_2[Lum_lne1_med_mc_2_indx] = np.nan
    Lum_lne2_avg_mc_1[Lum_lne2_avg_mc_1_indx] = np.nan
    Lum_lne2_med_mc_1[Lum_lne2_med_mc_1_indx] = np.nan
    Lum_lne2_avg_mc_2[Lum_lne2_avg_mc_2_indx] = np.nan
    Lum_lne2_med_mc_2[Lum_lne2_med_mc_2_indx] = np.nan

    Lum_lne12_avg_mc_1 = Lum_lne1_avg_mc_1 / Lum_lne2_avg_mc_1 
    Lum_lne12_avg_mc_2 = Lum_lne1_avg_mc_2 / Lum_lne2_avg_mc_2 
    Lum_lne12_med_mc_1 = Lum_lne1_med_mc_1 / Lum_lne2_med_mc_1 
    Lum_lne12_med_mc_2 = Lum_lne1_med_mc_2 / Lum_lne2_med_mc_2 

    ###################################################################1-SGM###################################################################
    z_lne1_avg_e_l_sgm1_mc,z_lne1_avg_e_s_sgm1_mc     = np.nanpercentile(z_z_lne1_mc,15.9),np.nanpercentile(z_z_lne1_mc,84.1)
    z_lne1_med_e_l_sgm1_mc,z_lne1_med_e_s_sgm1_mc     = np.nanpercentile(z_z_lne1_mc,15.9),np.nanpercentile(z_z_lne1_mc,84.1)
    z_lne2_avg_e_l_sgm1_mc,z_lne2_avg_e_s_sgm1_mc     = np.nanpercentile(z_z_lne2_mc,15.9),np.nanpercentile(z_z_lne2_mc,84.1)
    z_lne2_med_e_l_sgm1_mc,z_lne2_med_e_s_sgm1_mc     = np.nanpercentile(z_z_lne2_mc,15.9),np.nanpercentile(z_z_lne2_mc,84.1)

    v_lne1_avg_e_l_sgm1_mc,v_lne1_avg_e_s_sgm1_mc     = np.nanpercentile(v_z_lne1_mc,15.9),np.nanpercentile(v_z_lne1_mc,84.1)
    v_lne1_med_e_l_sgm1_mc,v_lne1_med_e_s_sgm1_mc     = np.nanpercentile(v_z_lne1_mc,15.9),np.nanpercentile(v_z_lne1_mc,84.1)
    v_lne2_avg_e_l_sgm1_mc,v_lne2_avg_e_s_sgm1_mc     = np.nanpercentile(v_z_lne2_mc,15.9),np.nanpercentile(v_z_lne2_mc,84.1)
    v_lne2_med_e_l_sgm1_mc,v_lne2_med_e_s_sgm1_mc     = np.nanpercentile(v_z_lne2_mc,15.9),np.nanpercentile(v_z_lne2_mc,84.1)

    s_lne1_avg_e_l_sgm1_mc,s_lne1_avg_e_s_sgm1_mc     = np.nanpercentile(s_lne1_avg_mc,15.9),np.nanpercentile(s_lne1_avg_mc,84.1)
    s_lne1_med_e_l_sgm1_mc,s_lne1_med_e_s_sgm1_mc     = np.nanpercentile(s_lne1_med_mc,15.9),np.nanpercentile(s_lne1_med_mc,84.1)
    s_lne2_avg_e_l_sgm1_mc,s_lne2_avg_e_s_sgm1_mc     = np.nanpercentile(s_lne2_avg_mc,15.9),np.nanpercentile(s_lne2_avg_mc,84.1)
    s_lne2_med_e_l_sgm1_mc,s_lne2_med_e_s_sgm1_mc     = np.nanpercentile(s_lne2_med_mc,15.9),np.nanpercentile(s_lne2_med_mc,84.1)
    
    Lum_lne1_avg_e_l_sgm1_mc_1,Lum_lne1_avg_e_s_sgm1_mc_1   = np.nanpercentile(Lum_lne1_avg_mc_1,15.9),np.nanpercentile(Lum_lne1_avg_mc_1,84.1)
    Lum_lne1_med_e_l_sgm1_mc_1,Lum_lne1_med_e_s_sgm1_mc_1   = np.nanpercentile(Lum_lne1_med_mc_1,15.9),np.nanpercentile(Lum_lne1_med_mc_1,84.1)
    Lum_lne1_avg_e_l_sgm1_mc_2,Lum_lne1_avg_e_s_sgm1_mc_2   = np.nanpercentile(Lum_lne1_avg_mc_2,15.9),np.nanpercentile(Lum_lne1_avg_mc_2,84.1)#log
    Lum_lne1_med_e_l_sgm1_mc_2,Lum_lne1_med_e_s_sgm1_mc_2   = np.nanpercentile(Lum_lne1_med_mc_2,15.9),np.nanpercentile(Lum_lne1_med_mc_2,84.1)#log
    Lum_lne2_avg_e_l_sgm1_mc_1,Lum_lne2_avg_e_s_sgm1_mc_1   = np.nanpercentile(Lum_lne2_avg_mc_1,15.9),np.nanpercentile(Lum_lne2_avg_mc_1,84.1)
    Lum_lne2_med_e_l_sgm1_mc_1,Lum_lne2_med_e_s_sgm1_mc_1   = np.nanpercentile(Lum_lne2_med_mc_1,15.9),np.nanpercentile(Lum_lne2_med_mc_1,84.1)
    Lum_lne2_avg_e_l_sgm1_mc_2,Lum_lne2_avg_e_s_sgm1_mc_2   = np.nanpercentile(Lum_lne2_avg_mc_2,15.9),np.nanpercentile(Lum_lne2_avg_mc_2,84.1)#log
    Lum_lne2_med_e_l_sgm1_mc_2,Lum_lne2_med_e_s_sgm1_mc_2   = np.nanpercentile(Lum_lne2_med_mc_2,15.9),np.nanpercentile(Lum_lne2_med_mc_2,84.1)#log

    z_lne12_avg_e_l_sgm1_mc,z_lne12_avg_e_s_sgm1_mc   = np.nanpercentile(z_lne12_avg_mc    ,15.9),np.nanpercentile(z_lne12_avg_mc    ,84.1)
    z_lne12_med_e_l_sgm1_mc,z_lne12_med_e_s_sgm1_mc   = np.nanpercentile(z_lne12_med_mc    ,15.9),np.nanpercentile(z_lne12_med_mc    ,84.1)
    v_lne12_avg_e_l_sgm1_mc,v_lne12_avg_e_s_sgm1_mc   = np.nanpercentile(v_lne12_avg_mc    ,15.9),np.nanpercentile(v_lne12_avg_mc    ,84.1)
    v_lne12_med_e_l_sgm1_mc,v_lne12_med_e_s_sgm1_mc   = np.nanpercentile(v_lne12_med_mc    ,15.9),np.nanpercentile(v_lne12_med_mc    ,84.1)
    s_lne12_avg_e_l_sgm1_mc,s_lne12_avg_e_s_sgm1_mc   = np.nanpercentile(s_lne12_avg_mc    ,15.9),np.nanpercentile(s_lne12_avg_mc    ,84.1)
    s_lne12_med_e_l_sgm1_mc,s_lne12_med_e_s_sgm1_mc   = np.nanpercentile(s_lne12_med_mc    ,15.9),np.nanpercentile(s_lne12_med_mc    ,84.1)

    Lum_lne12_avg_e_l_sgm1_mc_1,Lum_lne12_avg_e_s_sgm1_mc_1 = np.nanpercentile(Lum_lne12_avg_mc_1,15.9),np.nanpercentile(Lum_lne12_avg_mc_1,84.1)
    Lum_lne12_med_e_l_sgm1_mc_1,Lum_lne12_med_e_s_sgm1_mc_1 = np.nanpercentile(Lum_lne12_med_mc_1,15.9),np.nanpercentile(Lum_lne12_med_mc_1,84.1)
    Lum_lne12_avg_e_l_sgm1_mc_2,Lum_lne12_avg_e_s_sgm1_mc_2 = np.nanpercentile(Lum_lne12_avg_mc_2,15.9),np.nanpercentile(Lum_lne12_avg_mc_2,84.1)#log
    Lum_lne12_med_e_l_sgm1_mc_2,Lum_lne12_med_e_s_sgm1_mc_2 = np.nanpercentile(Lum_lne12_med_mc_2,15.9),np.nanpercentile(Lum_lne12_med_mc_2,84.1)#log
    ###################################################################1-SGM###################################################################
    ###################################################################2-SGM###################################################################

    z_lne1_avg_e_l_sgm2_mc,z_lne1_avg_e_s_sgm2_mc     = np.nanpercentile(z_z_lne1_mc,2.30),np.nanpercentile(z_z_lne1_mc,97.7)
    z_lne1_med_e_l_sgm2_mc,z_lne1_med_e_s_sgm2_mc     = np.nanpercentile(z_z_lne1_mc,2.30),np.nanpercentile(z_z_lne1_mc,97.7)
    z_lne2_avg_e_l_sgm2_mc,z_lne2_avg_e_s_sgm2_mc     = np.nanpercentile(z_z_lne2_mc,2.30),np.nanpercentile(z_z_lne2_mc,97.7)
    z_lne2_med_e_l_sgm2_mc,z_lne2_med_e_s_sgm2_mc     = np.nanpercentile(z_z_lne2_mc,2.30),np.nanpercentile(z_z_lne2_mc,97.7)

    v_lne1_avg_e_l_sgm2_mc,v_lne1_avg_e_s_sgm2_mc     = np.nanpercentile(v_z_lne1_mc,2.30),np.nanpercentile(v_z_lne1_mc,97.7)
    v_lne1_med_e_l_sgm2_mc,v_lne1_med_e_s_sgm2_mc     = np.nanpercentile(v_z_lne1_mc,2.30),np.nanpercentile(v_z_lne1_mc,97.7)
    v_lne2_avg_e_l_sgm2_mc,v_lne2_avg_e_s_sgm2_mc     = np.nanpercentile(v_z_lne2_mc,2.30),np.nanpercentile(v_z_lne2_mc,97.7)
    v_lne2_med_e_l_sgm2_mc,v_lne2_med_e_s_sgm2_mc     = np.nanpercentile(v_z_lne2_mc,2.30),np.nanpercentile(v_z_lne2_mc,97.7)

    s_lne1_avg_e_l_sgm2_mc,s_lne1_avg_e_s_sgm2_mc     = np.nanpercentile(s_lne1_avg_mc,2.30),np.nanpercentile(s_lne1_avg_mc,97.7)
    s_lne1_med_e_l_sgm2_mc,s_lne1_med_e_s_sgm2_mc     = np.nanpercentile(s_lne1_med_mc,2.30),np.nanpercentile(s_lne1_med_mc,97.7)
    s_lne2_avg_e_l_sgm2_mc,s_lne2_avg_e_s_sgm2_mc     = np.nanpercentile(s_lne2_avg_mc,2.30),np.nanpercentile(s_lne2_avg_mc,97.7)
    s_lne2_med_e_l_sgm2_mc,s_lne2_med_e_s_sgm2_mc     = np.nanpercentile(s_lne2_med_mc,2.30),np.nanpercentile(s_lne2_med_mc,97.7)
    #
    Lum_lne1_avg_e_l_sgm2_mc_1,Lum_lne1_avg_e_s_sgm2_mc_1   = np.nanpercentile(Lum_lne1_avg_mc_1,2.30),np.nanpercentile(Lum_lne1_avg_mc_1,97.7)
    Lum_lne1_med_e_l_sgm2_mc_1,Lum_lne1_med_e_s_sgm2_mc_1   = np.nanpercentile(Lum_lne1_med_mc_1,2.30),np.nanpercentile(Lum_lne1_med_mc_1,97.7)
    Lum_lne1_avg_e_l_sgm2_mc_2,Lum_lne1_avg_e_s_sgm2_mc_2   = np.nanpercentile(Lum_lne1_avg_mc_2,2.30),np.nanpercentile(Lum_lne1_avg_mc_2,97.7)#log
    Lum_lne1_med_e_l_sgm2_mc_2,Lum_lne1_med_e_s_sgm2_mc_2   = np.nanpercentile(Lum_lne1_med_mc_2,2.30),np.nanpercentile(Lum_lne1_med_mc_2,97.7)#log
    #
    Lum_lne2_avg_e_l_sgm2_mc_1,Lum_lne2_avg_e_s_sgm2_mc_1   = np.nanpercentile(Lum_lne2_avg_mc_1,2.30),np.nanpercentile(Lum_lne2_avg_mc_1,97.7)
    Lum_lne2_med_e_l_sgm2_mc_1,Lum_lne2_med_e_s_sgm2_mc_1   = np.nanpercentile(Lum_lne2_med_mc_1,2.30),np.nanpercentile(Lum_lne2_med_mc_1,97.7)
    Lum_lne2_avg_e_l_sgm2_mc_2,Lum_lne2_avg_e_s_sgm2_mc_2   = np.nanpercentile(Lum_lne2_avg_mc_2,2.30),np.nanpercentile(Lum_lne2_avg_mc_2,97.7)#log
    Lum_lne2_med_e_l_sgm2_mc_2,Lum_lne2_med_e_s_sgm2_mc_2   = np.nanpercentile(Lum_lne2_med_mc_2,2.30),np.nanpercentile(Lum_lne2_med_mc_2,97.7)#log
    #
    z_lne12_avg_e_l_sgm2_mc,z_lne12_avg_e_s_sgm2_mc   = np.nanpercentile(z_lne12_avg_mc    ,2.30),np.nanpercentile(z_lne12_avg_mc    ,97.7)
    z_lne12_med_e_l_sgm2_mc,z_lne12_med_e_s_sgm2_mc   = np.nanpercentile(z_lne12_med_mc    ,2.30),np.nanpercentile(z_lne12_med_mc    ,97.7)
    v_lne12_avg_e_l_sgm2_mc,v_lne12_avg_e_s_sgm2_mc   = np.nanpercentile(v_lne12_avg_mc    ,2.30),np.nanpercentile(v_lne12_avg_mc    ,97.7)
    v_lne12_med_e_l_sgm2_mc,v_lne12_med_e_s_sgm2_mc   = np.nanpercentile(v_lne12_med_mc    ,2.30),np.nanpercentile(v_lne12_med_mc    ,97.7)
    s_lne12_avg_e_l_sgm2_mc,s_lne12_avg_e_s_sgm2_mc   = np.nanpercentile(s_lne12_avg_mc    ,2.30),np.nanpercentile(s_lne12_avg_mc    ,97.7)
    s_lne12_med_e_l_sgm2_mc,s_lne12_med_e_s_sgm2_mc   = np.nanpercentile(s_lne12_med_mc    ,2.30),np.nanpercentile(s_lne12_med_mc    ,97.7)

    Lum_lne12_avg_e_l_sgm2_mc_1,Lum_lne12_avg_e_s_sgm2_mc_1 = np.nanpercentile(Lum_lne12_avg_mc_1,2.30),np.nanpercentile(Lum_lne12_avg_mc_1,97.7)
    Lum_lne12_med_e_l_sgm2_mc_1,Lum_lne12_med_e_s_sgm2_mc_1 = np.nanpercentile(Lum_lne12_med_mc_1,2.30),np.nanpercentile(Lum_lne12_med_mc_1,97.7)
    Lum_lne12_avg_e_l_sgm2_mc_2,Lum_lne12_avg_e_s_sgm2_mc_2 = np.nanpercentile(Lum_lne12_avg_mc_2,2.30),np.nanpercentile(Lum_lne12_avg_mc_2,97.7)#log
    Lum_lne12_med_e_l_sgm2_mc_2,Lum_lne12_med_e_s_sgm2_mc_2 = np.nanpercentile(Lum_lne12_med_mc_2,2.30),np.nanpercentile(Lum_lne12_med_mc_2,97.7)#log
    ###################################################################2-SGM###################################################################

    ###################################################################3-SGM###################################################################
    z_lne1_avg_e_l_sgm3_mc,z_lne1_avg_e_s_sgm3_mc     = np.nanpercentile(z_z_lne1_mc,0.20),np.nanpercentile(z_z_lne1_mc,99.8)
    z_lne1_med_e_l_sgm3_mc,z_lne1_med_e_s_sgm3_mc     = np.nanpercentile(z_z_lne1_mc,0.20),np.nanpercentile(z_z_lne1_mc,99.8)
    z_lne2_avg_e_l_sgm3_mc,z_lne2_avg_e_s_sgm3_mc     = np.nanpercentile(z_z_lne2_mc,0.20),np.nanpercentile(z_z_lne2_mc,99.8)
    z_lne2_med_e_l_sgm3_mc,z_lne2_med_e_s_sgm3_mc     = np.nanpercentile(z_z_lne2_mc,0.20),np.nanpercentile(z_z_lne2_mc,99.8)

    v_lne1_avg_e_l_sgm3_mc,v_lne1_avg_e_s_sgm3_mc     = np.nanpercentile(v_z_lne1_mc,0.20),np.nanpercentile(v_z_lne1_mc,99.8)
    v_lne1_med_e_l_sgm3_mc,v_lne1_med_e_s_sgm3_mc     = np.nanpercentile(v_z_lne1_mc,0.20),np.nanpercentile(v_z_lne1_mc,99.8)
    v_lne2_avg_e_l_sgm3_mc,v_lne2_avg_e_s_sgm3_mc     = np.nanpercentile(v_z_lne2_mc,0.20),np.nanpercentile(v_z_lne2_mc,99.8)
    v_lne2_med_e_l_sgm3_mc,v_lne2_med_e_s_sgm3_mc     = np.nanpercentile(v_z_lne2_mc,0.20),np.nanpercentile(v_z_lne2_mc,99.8)

    s_lne1_avg_e_l_sgm3_mc,s_lne1_avg_e_s_sgm3_mc     = np.nanpercentile(s_lne1_avg_mc,0.20),np.nanpercentile(s_lne1_avg_mc,99.8)
    s_lne1_med_e_l_sgm3_mc,s_lne1_med_e_s_sgm3_mc     = np.nanpercentile(s_lne1_med_mc,0.20),np.nanpercentile(s_lne1_med_mc,99.8)
    s_lne2_avg_e_l_sgm3_mc,s_lne2_avg_e_s_sgm3_mc     = np.nanpercentile(s_lne2_avg_mc,0.20),np.nanpercentile(s_lne2_avg_mc,99.8)
    s_lne2_med_e_l_sgm3_mc,s_lne2_med_e_s_sgm3_mc     = np.nanpercentile(s_lne2_med_mc,0.20),np.nanpercentile(s_lne2_med_mc,99.8)
    #
    Lum_lne1_avg_e_l_sgm3_mc_1,Lum_lne1_avg_e_s_sgm3_mc_1   = np.nanpercentile(Lum_lne1_avg_mc_1,0.20),np.nanpercentile(Lum_lne1_avg_mc_1,99.8)
    Lum_lne1_med_e_l_sgm3_mc_1,Lum_lne1_med_e_s_sgm3_mc_1   = np.nanpercentile(Lum_lne1_med_mc_1,0.20),np.nanpercentile(Lum_lne1_med_mc_1,99.8)
    Lum_lne1_avg_e_l_sgm3_mc_2,Lum_lne1_avg_e_s_sgm3_mc_2   = np.nanpercentile(Lum_lne1_avg_mc_2,0.20),np.nanpercentile(Lum_lne1_avg_mc_2,99.8)#log
    Lum_lne1_med_e_l_sgm3_mc_2,Lum_lne1_med_e_s_sgm3_mc_2   = np.nanpercentile(Lum_lne1_med_mc_2,0.20),np.nanpercentile(Lum_lne1_med_mc_2,99.8)#log
    
    Lum_lne2_avg_e_l_sgm3_mc_1,Lum_lne2_avg_e_s_sgm3_mc_1   = np.nanpercentile(Lum_lne2_avg_mc_1,0.20),np.nanpercentile(Lum_lne2_avg_mc_1,99.8)
    Lum_lne2_med_e_l_sgm3_mc_1,Lum_lne2_med_e_s_sgm3_mc_1   = np.nanpercentile(Lum_lne2_med_mc_1,0.20),np.nanpercentile(Lum_lne2_med_mc_1,99.8)
    Lum_lne2_avg_e_l_sgm3_mc_2,Lum_lne2_avg_e_s_sgm3_mc_2   = np.nanpercentile(Lum_lne2_avg_mc_2,0.20),np.nanpercentile(Lum_lne2_avg_mc_2,99.8)#log
    Lum_lne2_med_e_l_sgm3_mc_2,Lum_lne2_med_e_s_sgm3_mc_2   = np.nanpercentile(Lum_lne2_med_mc_2,0.20),np.nanpercentile(Lum_lne2_med_mc_2,99.8)#log

    z_lne12_avg_e_l_sgm3_mc ,z_lne12_avg_e_s_sgm3_mc  = np.nanpercentile(z_lne12_avg_mc,0.20),np.nanpercentile(z_lne12_avg_mc    ,99.8)
    z_lne12_med_e_l_sgm3_mc ,z_lne12_med_e_s_sgm3_mc  = np.nanpercentile(z_lne12_med_mc,0.20),np.nanpercentile(z_lne12_med_mc    ,99.8)
    v_lne12_avg_e_l_sgm3_mc ,v_lne12_avg_e_s_sgm3_mc  = np.nanpercentile(v_lne12_avg_mc,0.20),np.nanpercentile(v_lne12_avg_mc    ,99.8)
    v_lne12_med_e_l_sgm3_mc ,v_lne12_med_e_s_sgm3_mc  = np.nanpercentile(v_lne12_med_mc,0.20),np.nanpercentile(v_lne12_med_mc    ,99.8)
    s_lne12_avg_e_l_sgm3_mc ,s_lne12_avg_e_s_sgm3_mc  = np.nanpercentile(s_lne12_avg_mc,0.20),np.nanpercentile(s_lne12_avg_mc    ,99.8)
    s_lne12_med_e_l_sgm3_mc ,s_lne12_med_e_s_sgm3_mc  = np.nanpercentile(s_lne12_med_mc,0.20),np.nanpercentile(s_lne12_med_mc    ,99.8)

    Lum_lne12_avg_e_l_sgm3_mc_1,Lum_lne12_avg_e_s_sgm3_mc_1 = np.nanpercentile(Lum_lne12_avg_mc_1,0.20),np.nanpercentile(Lum_lne12_avg_mc_1,99.8)
    Lum_lne12_med_e_l_sgm3_mc_1,Lum_lne12_med_e_s_sgm3_mc_1 = np.nanpercentile(Lum_lne12_med_mc_1,0.20),np.nanpercentile(Lum_lne12_med_mc_1,99.8)
    Lum_lne12_avg_e_l_sgm3_mc_2,Lum_lne12_avg_e_s_sgm3_mc_2 = np.nanpercentile(Lum_lne12_avg_mc_2,0.20),np.nanpercentile(Lum_lne12_avg_mc_2,99.8)#log
    Lum_lne12_med_e_l_sgm3_mc_2,Lum_lne12_med_e_s_sgm3_mc_2 = np.nanpercentile(Lum_lne12_med_mc_2,0.20),np.nanpercentile(Lum_lne12_med_mc_2,99.8)#log
    ###################################################################3-SGM###################################################################
    #######################################################ADD-MC-ERRORS-TO-FITS-HEADER########################################################

    print
    print ('1: ',spec_file_plt_lne1_avg_mcmc)
    print ('1: ',spec_file_plt_lne1_med_mcmc)
    print ('2: ',spec_file_plt_lne2_avg_mcmc)
    print ('2: ',spec_file_plt_lne2_med_mcmc)
    print

    HEADER_T = [
            'MSR_TYP' + str(method),'MCR_NMB' + str(method),
            'ZMCE_1L' + str(method),'ZMCE_1H' + str(method),'ZMCE_2L' + str(method),'ZMCE_2H' + str(method),'ZMCE_3L' + str(method),'ZMCE_3H'+str(method),
            'VMCE_1L' + str(method),'VMCE_1H' + str(method),'VMCE_2L' + str(method),'VMCE_2H' + str(method),'VMCE_3L' + str(method),'VMCE_3H'+str(method),
            'FMCE_1L' + str(method),'FMCE_1H' + str(method),'FMCE_2L' + str(method),'FMCE_2H' + str(method),'FMCE_3L' + str(method),'FMCE_3H'+str(method),
            'LMCE_1L' + str(method),'LMCE_1H' + str(method),'LMCE_2L' + str(method),'LMCE_2H' + str(method),'LMCE_3L' + str(method),'LMCE_3H'+str(method),
            'LLME_1L' + str(method),'LLME_1H' + str(method),'LLME_2L' + str(method),'LLME_2H' + str(method),'LLME_3L' + str(method),'LLME_3H'+str(method),
            'ZMER_1L' + str(method),'ZMER_1H' + str(method),'ZMER_2L' + str(method),'ZMER_2H' + str(method),'ZMER_3L' + str(method),'ZMER_3H'+str(method),
            'VMER_1L' + str(method),'VMER_1H' + str(method),'VMER_2L' + str(method),'VMER_2H' + str(method),'VMER_3L' + str(method),'VMER_3H'+str(method),
            'FMER_1L' + str(method),'FMER_1H' + str(method),'FMER_2L' + str(method),'FMER_2H' + str(method),'FMER_3L' + str(method),'FMER_3H'+str(method),
            'LMER_1L' + str(method),'LMER_1H' + str(method),'LMER_2L' + str(method),'LMER_2H' + str(method),'LMER_3L' + str(method),'LMER_3H'+str(method),
            'LLMR_1L' + str(method),'LLMR_1H' + str(method),'LLMR_2L' + str(method),'LLMR_2H' + str(method),'LLMR_3L' + str(method),'LLMR_3H'+str(method) 
            ]

    HEADER_C = ['Tot Flx Msr Type: 1-sum 2-1GF 3-2DGF 4-2GF-fwhm','MCE Repetitions: '+ str(iterations_mc),
            'MCE Redshift 1 sgm lw lmt 15.9 pct','MCE Redshift 1 sgm hg lmt 84.1 pct',
            'MCE Redshift 2 sgm lw lmt 2.30 pct','MCE Redshift 2 sgm hg lmt 97.7 pct',
            'MCE Redshift 3 sgm lw lmt 0.20 pct','MCE Redshift 3 sgm hg lmt 99.8 pct',
            'MCE ' +str(sbsms) + ' 1 sgm lw lmt 15.9 pct','MCE ' +str(sbsms) + ' 1 sgm hg lmt 84.1 pct',
            'MCE ' +str(sbsms) + ' 2 sgm lw lmt 2.30 pct','MCE ' +str(sbsms) + ' 2 sgm hg lmt 97.7 pct',
            'MCE ' +str(sbsms) + ' 3 sgm lw lmt 0.20 pct','MCE ' +str(sbsms) + ' 3 sgm hg lmt 99.8 pct',
            'MCE Flx 1 sgm lw lmt 15.9 pct','MCE Flx 1 sgm hg lmt 84.1 pct',
            'MCE Flx 2 sgm lw lmt 2.30 pct','MCE Flx 2 sgm hg lmt 97.7 pct',
            'MCE Flx 3 sgm lw lmt 0.20 pct','MCE Flx 3 sgm hg lmt 99.8 pct',
            'MCE Lum 1 sgm lw lmt 15.9 pct','MCE Lum 1 sgm hg lmt 84.1 pct',
            'MCE Lum 2 sgm lw lmt 2.30 pct','MCE Lum 2 sgm hg lmt 97.7 pct',
            'MCE Lum 3 sgm lw lmt 0.20 pct','MCE Lum 3 sgm hg lmt 99.8 pct',
            'MCE Lg Lum 1 sgm lw lmt 15.9 pct','MCE Lg Lum 1 sgm hg lmt 84.1 pct',
            'MCE Lg Lum 2 sgm lw lmt 2.30 pct','MCE Lg Lum 2 sgm hg lmt 97.7 pct',
            'MCE Lg Lum 3 sgm lw lmt 0.20 pct','MCE Lg Lum 3 sgm hg lmt 99.8 pct',
            'MCE Redshift 12R 1 sgm lw lmt 15.9 pct','MCE Redshift 12R 1 sgm hg lmt 84.1 pct',
            'MCE Redshift 12R 2 sgm lw lmt 2.30 pct','MCE Redshift 12R 2 sgm hg lmt 97.7 pct',
            'MCE Redshift 12R 3 sgm lw lmt 0.20 pct','MCE Redshift 12R 3 sgm hg lmt 99.8 pct',
            'MCE ' + str(sbsms) + ' 12R 1 sgm lw lmt 15.9 pct','MCE ' + str(sbsms) + ' 12R 1 sgm hg lmt 84.1 pct',
            'MCE ' + str(sbsms) + ' 12R 2 sgm lw lmt 2.30 pct','MCE ' + str(sbsms) + ' 12R 2 sgm hg lmt 97.7 pct',
            'MCE ' + str(sbsms) + ' 12R 3 sgm lw lmt 0.20 pct','MCE ' + str(sbsms) + ' 12R 3 sgm hg lmt 99.8 pct',
            'MCE Flux 12R 1 sgm lw lmt 15.9 pct','MCE Flux 12R 1 sgm hg lmt 84.1 pct',
            'MCE Flux 12R 2 sgm lw lmt 2.30 pct','MCE Flux 12R 2 sgm hg lmt 97.7 pct',
            'MCE Flux 12R 3 sgm lw lmt 0.20 pct','MCE Flux 12R 3 sgm hg lmt 99.8 pct',
            'MCE Lum 12R 1 sgm lw lmt 15.9 pct','MCE Lum 12R 1 sgm hg lmt 84.1 pct',
            'MCE Lum 12R 2 sgm lw lmt 2.30 pct','MCE Lum 12R 2 sgm hg lmt 97.7 pct',
            'MCE Lum 12R 3 sgm lw lmt 0.20 pct','MCE Lum 12R 3 sgm hg lmt 99.8 pct',
            'MCE Lg Lum 12R 1 sgm lw lmt 15.9 pct','MCE Lg Lum 12R 1 sgm hg lmt 84.1 pct',
            'MCE Lg Lum 12R 2 sgm lw lmt 2.30 pct','MCE Lg Lum 12R 2 sgm hg lmt 97.7 pct',
            'MCE Lg Lum 12R 3 sgm lw lmt 0.20 pct','MCE Lg Lum 12R 3 sgm hg lmt 99.8 pct'
            ]
            
    HEADER_V = [method,iterations_mc,
            z_lne1_avg_e_l_sgm1_mc,z_lne1_avg_e_s_sgm1_mc,z_lne1_avg_e_l_sgm2_mc,z_lne1_avg_e_s_sgm2_mc,z_lne1_avg_e_l_sgm3_mc,z_lne1_avg_e_s_sgm3_mc,
            v_lne1_avg_e_l_sgm1_mc,v_lne1_avg_e_s_sgm1_mc,v_lne1_avg_e_l_sgm2_mc,v_lne1_avg_e_s_sgm2_mc,v_lne1_avg_e_l_sgm3_mc,v_lne1_avg_e_s_sgm3_mc,
            s_lne1_avg_e_l_sgm1_mc,s_lne1_avg_e_s_sgm1_mc,s_lne1_avg_e_l_sgm2_mc,s_lne1_avg_e_s_sgm2_mc,s_lne1_avg_e_l_sgm3_mc,s_lne1_avg_e_s_sgm3_mc,
            Lum_lne1_avg_e_l_sgm1_mc_1,Lum_lne1_avg_e_s_sgm1_mc_1,Lum_lne1_avg_e_l_sgm2_mc_1,Lum_lne1_avg_e_s_sgm2_mc_1,Lum_lne1_avg_e_l_sgm3_mc_1,Lum_lne1_avg_e_s_sgm3_mc_1,
            Lum_lne1_avg_e_l_sgm1_mc_2,Lum_lne1_avg_e_s_sgm1_mc_2,Lum_lne1_avg_e_l_sgm2_mc_2,Lum_lne1_avg_e_s_sgm2_mc_2,Lum_lne1_avg_e_l_sgm3_mc_2,Lum_lne1_avg_e_s_sgm3_mc_2,
            z_lne12_avg_e_l_sgm1_mc,z_lne12_avg_e_s_sgm1_mc,z_lne12_avg_e_l_sgm2_mc,z_lne12_avg_e_s_sgm2_mc,z_lne12_avg_e_l_sgm3_mc,z_lne12_avg_e_s_sgm3_mc,
            v_lne12_avg_e_l_sgm1_mc,v_lne12_avg_e_s_sgm1_mc,v_lne12_avg_e_l_sgm2_mc,v_lne12_avg_e_s_sgm2_mc,v_lne12_avg_e_l_sgm3_mc,v_lne12_avg_e_s_sgm3_mc,
            s_lne12_avg_e_l_sgm1_mc,s_lne12_avg_e_s_sgm1_mc,s_lne12_avg_e_l_sgm2_mc,s_lne12_avg_e_s_sgm2_mc,s_lne12_avg_e_l_sgm3_mc,s_lne12_avg_e_s_sgm3_mc,
            Lum_lne12_avg_e_l_sgm1_mc_1,Lum_lne12_avg_e_s_sgm1_mc_1,Lum_lne12_avg_e_l_sgm2_mc_1,Lum_lne12_avg_e_s_sgm2_mc_1,Lum_lne12_avg_e_l_sgm3_mc_1,Lum_lne12_avg_e_s_sgm3_mc_1,
            Lum_lne12_avg_e_l_sgm1_mc_2,Lum_lne12_avg_e_s_sgm1_mc_2,Lum_lne12_avg_e_l_sgm2_mc_2,Lum_lne12_avg_e_s_sgm2_mc_2,Lum_lne12_avg_e_l_sgm3_mc_2,Lum_lne12_avg_e_s_sgm3_mc_2
            ]

    [Header_Add(spec_file_plt_lne1_avg_mcmc,head_nme,head_val,header_comment=head_com) for head_nme,head_val,head_com in zip(HEADER_T,HEADER_V,HEADER_C)]

    HEADER_V = [method,iterations_mc,
            z_lne1_med_e_l_sgm1_mc,z_lne1_med_e_s_sgm1_mc,z_lne1_med_e_l_sgm2_mc,z_lne1_med_e_s_sgm2_mc,z_lne1_med_e_l_sgm3_mc,z_lne1_med_e_s_sgm3_mc,
            v_lne1_med_e_l_sgm1_mc,v_lne1_med_e_s_sgm1_mc,v_lne1_med_e_l_sgm2_mc,v_lne1_med_e_s_sgm2_mc,v_lne1_med_e_l_sgm3_mc,v_lne1_med_e_s_sgm3_mc,
            s_lne1_med_e_l_sgm1_mc,s_lne1_med_e_s_sgm1_mc,s_lne1_med_e_l_sgm2_mc,s_lne1_med_e_s_sgm2_mc,s_lne1_med_e_l_sgm3_mc,s_lne1_med_e_s_sgm3_mc,
            Lum_lne1_med_e_l_sgm1_mc_1,Lum_lne1_med_e_s_sgm1_mc_1,Lum_lne1_med_e_l_sgm2_mc_1,Lum_lne1_med_e_s_sgm2_mc_1,Lum_lne1_med_e_l_sgm3_mc_1,Lum_lne1_med_e_s_sgm3_mc_1,
            Lum_lne1_med_e_l_sgm1_mc_2,Lum_lne1_med_e_s_sgm1_mc_2,Lum_lne1_med_e_l_sgm2_mc_2,Lum_lne1_med_e_s_sgm2_mc_2,Lum_lne1_med_e_l_sgm3_mc_2,Lum_lne1_med_e_s_sgm3_mc_2,
            z_lne12_med_e_l_sgm1_mc,z_lne12_med_e_s_sgm1_mc,z_lne12_med_e_l_sgm2_mc,z_lne12_med_e_s_sgm2_mc,z_lne12_med_e_l_sgm3_mc,z_lne12_med_e_s_sgm3_mc,
            v_lne12_med_e_l_sgm1_mc,v_lne12_med_e_s_sgm1_mc,v_lne12_med_e_l_sgm2_mc,v_lne12_med_e_s_sgm2_mc,v_lne12_med_e_l_sgm3_mc,v_lne12_med_e_s_sgm3_mc,
            s_lne12_med_e_l_sgm1_mc,s_lne12_med_e_s_sgm1_mc,s_lne12_med_e_l_sgm2_mc,s_lne12_med_e_s_sgm2_mc,s_lne12_med_e_l_sgm3_mc,s_lne12_med_e_s_sgm3_mc,
            Lum_lne12_med_e_l_sgm1_mc_1,Lum_lne12_med_e_s_sgm1_mc_1,Lum_lne12_med_e_l_sgm2_mc_1,Lum_lne12_med_e_s_sgm2_mc_1,Lum_lne12_med_e_l_sgm3_mc_1,Lum_lne12_med_e_s_sgm3_mc_1,
            Lum_lne12_med_e_l_sgm1_mc_2,Lum_lne12_med_e_s_sgm1_mc_2,Lum_lne12_med_e_l_sgm2_mc_2,Lum_lne12_med_e_s_sgm2_mc_2,Lum_lne12_med_e_l_sgm3_mc_2,Lum_lne12_med_e_s_sgm3_mc_2
            ]
    [Header_Add(spec_file_plt_lne1_med_mcmc,head_nme,head_val,header_comment=head_com) for head_nme,head_val,head_com in zip(HEADER_T,HEADER_V,HEADER_C)]

    HEADER_V = [method,iterations_mc,
            z_lne2_avg_e_l_sgm1_mc,z_lne2_avg_e_s_sgm1_mc,z_lne2_avg_e_l_sgm2_mc,z_lne2_avg_e_s_sgm2_mc,z_lne2_avg_e_l_sgm3_mc,z_lne2_avg_e_s_sgm3_mc,
            v_lne2_avg_e_l_sgm1_mc,v_lne2_avg_e_s_sgm1_mc,v_lne2_avg_e_l_sgm2_mc,v_lne2_avg_e_s_sgm2_mc,v_lne2_avg_e_l_sgm3_mc,v_lne2_avg_e_s_sgm3_mc,
            s_lne2_avg_e_l_sgm1_mc,s_lne2_avg_e_s_sgm1_mc,s_lne2_avg_e_l_sgm2_mc,s_lne2_avg_e_s_sgm2_mc,s_lne2_avg_e_l_sgm3_mc,s_lne2_avg_e_s_sgm3_mc,
            Lum_lne2_avg_e_l_sgm1_mc_1,Lum_lne2_avg_e_s_sgm1_mc_1,Lum_lne2_avg_e_l_sgm2_mc_1,Lum_lne2_avg_e_s_sgm2_mc_1,Lum_lne2_avg_e_l_sgm3_mc_1,Lum_lne2_avg_e_s_sgm3_mc_1,
            Lum_lne2_avg_e_l_sgm1_mc_2,Lum_lne2_avg_e_s_sgm1_mc_2,Lum_lne2_avg_e_l_sgm2_mc_2,Lum_lne2_avg_e_s_sgm2_mc_2,Lum_lne2_avg_e_l_sgm3_mc_2,Lum_lne2_avg_e_s_sgm3_mc_2,
            z_lne12_avg_e_l_sgm1_mc,z_lne12_avg_e_s_sgm1_mc,z_lne12_avg_e_l_sgm2_mc,z_lne12_avg_e_s_sgm2_mc,z_lne12_avg_e_l_sgm3_mc,z_lne12_avg_e_s_sgm3_mc,
            v_lne12_avg_e_l_sgm1_mc,v_lne12_avg_e_s_sgm1_mc,v_lne12_avg_e_l_sgm2_mc,v_lne12_avg_e_s_sgm2_mc,v_lne12_avg_e_l_sgm3_mc,v_lne12_avg_e_s_sgm3_mc,
            s_lne12_avg_e_l_sgm1_mc,s_lne12_avg_e_s_sgm1_mc,s_lne12_avg_e_l_sgm2_mc,s_lne12_avg_e_s_sgm2_mc,s_lne12_avg_e_l_sgm3_mc,s_lne12_avg_e_s_sgm3_mc,
            Lum_lne12_avg_e_l_sgm1_mc_1,Lum_lne12_avg_e_s_sgm1_mc_1,Lum_lne12_avg_e_l_sgm2_mc_1,Lum_lne12_avg_e_s_sgm2_mc_1,Lum_lne12_avg_e_l_sgm3_mc_1,Lum_lne12_avg_e_s_sgm3_mc_1,
            Lum_lne12_avg_e_l_sgm1_mc_2,Lum_lne12_avg_e_s_sgm1_mc_2,Lum_lne12_avg_e_l_sgm2_mc_2,Lum_lne12_avg_e_s_sgm2_mc_2,Lum_lne12_avg_e_l_sgm3_mc_2,Lum_lne12_avg_e_s_sgm3_mc_2
            ]
    [Header_Add(spec_file_plt_lne2_avg_mcmc,head_nme,head_val,header_comment=head_com) for head_nme,head_val,head_com in zip(HEADER_T,HEADER_V,HEADER_C)]

    HEADER_V = [method,iterations_mc,
            z_lne2_med_e_l_sgm1_mc,z_lne2_med_e_s_sgm1_mc,z_lne2_med_e_l_sgm2_mc,z_lne2_med_e_s_sgm2_mc,z_lne2_med_e_l_sgm3_mc,z_lne2_med_e_s_sgm3_mc,
            v_lne2_med_e_l_sgm1_mc,v_lne2_med_e_s_sgm1_mc,v_lne2_med_e_l_sgm2_mc,v_lne2_med_e_s_sgm2_mc,v_lne2_med_e_l_sgm3_mc,v_lne2_med_e_s_sgm3_mc,
            s_lne2_med_e_l_sgm1_mc,s_lne2_med_e_s_sgm1_mc,s_lne2_med_e_l_sgm2_mc,s_lne2_med_e_s_sgm2_mc,s_lne2_med_e_l_sgm3_mc,s_lne2_med_e_s_sgm3_mc,
            Lum_lne2_med_e_l_sgm1_mc_1,Lum_lne2_med_e_s_sgm1_mc_1,Lum_lne2_med_e_l_sgm2_mc_1,Lum_lne2_med_e_s_sgm2_mc_1,Lum_lne2_med_e_l_sgm3_mc_1,Lum_lne2_med_e_s_sgm3_mc_1,
            Lum_lne2_med_e_l_sgm1_mc_2,Lum_lne2_med_e_s_sgm1_mc_2,Lum_lne2_med_e_l_sgm2_mc_2,Lum_lne2_med_e_s_sgm2_mc_2,Lum_lne2_med_e_l_sgm3_mc_2,Lum_lne2_med_e_s_sgm3_mc_2,
            z_lne12_med_e_l_sgm1_mc,z_lne12_med_e_s_sgm1_mc,z_lne12_med_e_l_sgm2_mc,z_lne12_med_e_s_sgm2_mc,z_lne12_med_e_l_sgm3_mc,z_lne12_med_e_s_sgm3_mc,
            v_lne12_med_e_l_sgm1_mc,v_lne12_med_e_s_sgm1_mc,v_lne12_med_e_l_sgm2_mc,v_lne12_med_e_s_sgm2_mc,v_lne12_med_e_l_sgm3_mc,v_lne12_med_e_s_sgm3_mc,
            s_lne12_med_e_l_sgm1_mc,s_lne12_med_e_s_sgm1_mc,s_lne12_med_e_l_sgm2_mc,s_lne12_med_e_s_sgm2_mc,s_lne12_med_e_l_sgm3_mc,s_lne12_med_e_s_sgm3_mc,
            Lum_lne12_med_e_l_sgm1_mc_1,Lum_lne12_med_e_s_sgm1_mc_1,Lum_lne12_med_e_l_sgm2_mc_1,Lum_lne12_med_e_s_sgm2_mc_1,Lum_lne12_med_e_l_sgm3_mc_1,Lum_lne12_med_e_s_sgm3_mc_1,
            Lum_lne12_med_e_l_sgm1_mc_2,Lum_lne12_med_e_s_sgm1_mc_2,Lum_lne12_med_e_l_sgm2_mc_2,Lum_lne12_med_e_s_sgm2_mc_2,Lum_lne12_med_e_l_sgm3_mc_2,Lum_lne12_med_e_s_sgm3_mc_2
            ]
    [Header_Add(spec_file_plt_lne2_med_mcmc,head_nme,head_val,header_comment=head_com) for head_nme,head_val,head_com in zip(HEADER_T,HEADER_V,HEADER_C)]

    ########################################################ADD-MC-ERRORS-TO-FITS-HEADER########################################################

    ###################################################################3-SGM###################################################################

    ################PLOT-ARRAYS######################
    Z1_AVG_E1_MC.append(z_lne1_avg_e_l_sgm1_mc)
    Z1_AVG_E2_MC.append(z_lne1_avg_e_s_sgm1_mc)
    Z2_AVG_E1_MC.append(z_lne2_avg_e_l_sgm1_mc)
    Z2_AVG_E2_MC.append(z_lne2_avg_e_s_sgm1_mc)

    Z1_MED_E1_MC.append(z_lne1_med_e_l_sgm1_mc)
    Z1_MED_E2_MC.append(z_lne1_med_e_s_sgm1_mc)
    Z2_MED_E1_MC.append(z_lne2_med_e_l_sgm1_mc)
    Z2_MED_E2_MC.append(z_lne2_med_e_s_sgm1_mc)

    V1_AVG_E1_MC.append(v_lne1_avg_e_l_sgm1_mc)
    V1_AVG_E2_MC.append(v_lne1_avg_e_s_sgm1_mc)
    V2_AVG_E1_MC.append(v_lne2_avg_e_l_sgm1_mc)
    V2_AVG_E2_MC.append(v_lne2_avg_e_s_sgm1_mc)

    V1_MED_E1_MC.append(v_lne1_med_e_l_sgm1_mc)
    V1_MED_E2_MC.append(v_lne1_med_e_s_sgm1_mc)
    V2_MED_E1_MC.append(v_lne2_med_e_l_sgm1_mc)
    V2_MED_E2_MC.append(v_lne2_med_e_s_sgm1_mc)

    S1_AVG_E1_MC.append(s_lne1_avg_e_l_sgm1_mc)
    S1_AVG_E2_MC.append(s_lne1_avg_e_s_sgm1_mc)
    S2_AVG_E1_MC.append(s_lne2_avg_e_l_sgm1_mc)
    S2_AVG_E2_MC.append(s_lne2_avg_e_s_sgm1_mc)

    S1_MED_E1_MC.append(s_lne1_med_e_l_sgm1_mc)
    S1_MED_E2_MC.append(s_lne1_med_e_s_sgm1_mc)
    S2_MED_E1_MC.append(s_lne2_med_e_l_sgm1_mc)
    S2_MED_E2_MC.append(s_lne2_med_e_s_sgm1_mc)

    L1_AVG_E1_MC_1.append(Lum_lne1_avg_e_l_sgm1_mc_1)
    L1_AVG_E2_MC_1.append(Lum_lne1_avg_e_s_sgm1_mc_1)
    L1_MED_E1_MC_1.append(Lum_lne1_med_e_l_sgm1_mc_1)
    L1_MED_E2_MC_1.append(Lum_lne1_med_e_s_sgm1_mc_1)

    L1_AVG_E1_MC_2.append(Lum_lne1_avg_e_l_sgm1_mc_2)
    L1_AVG_E2_MC_2.append(Lum_lne1_avg_e_s_sgm1_mc_2)
    L1_MED_E1_MC_2.append(Lum_lne1_med_e_l_sgm1_mc_2)
    L1_MED_E2_MC_2.append(Lum_lne1_med_e_s_sgm1_mc_2)

    L2_AVG_E1_MC_1.append(Lum_lne2_avg_e_l_sgm1_mc_1)
    L2_AVG_E2_MC_1.append(Lum_lne2_avg_e_s_sgm1_mc_1)
    L2_MED_E1_MC_1.append(Lum_lne2_med_e_l_sgm1_mc_1)
    L2_MED_E2_MC_1.append(Lum_lne2_med_e_s_sgm1_mc_1)

    L2_AVG_E1_MC_2.append(Lum_lne2_avg_e_l_sgm1_mc_2)
    L2_AVG_E2_MC_2.append(Lum_lne2_avg_e_s_sgm1_mc_2)
    L2_MED_E1_MC_2.append(Lum_lne2_med_e_l_sgm1_mc_2)
    L2_MED_E2_MC_2.append(Lum_lne2_med_e_s_sgm1_mc_2)

    S12_AVG_E1_MC.append(s_lne12_avg_e_l_sgm1_mc)
    S12_AVG_E2_MC.append(s_lne12_avg_e_s_sgm1_mc)
    S12_MED_E1_MC.append(s_lne12_med_e_l_sgm1_mc)
    S12_MED_E2_MC.append(s_lne12_med_e_s_sgm1_mc)

    L12_AVG_E1_MC_1.append(Lum_lne12_avg_e_l_sgm1_mc_1)
    L12_AVG_E2_MC_1.append(Lum_lne12_avg_e_s_sgm1_mc_1)
    L12_MED_E1_MC_1.append(Lum_lne12_med_e_l_sgm1_mc_1)
    L12_MED_E2_MC_1.append(Lum_lne12_med_e_s_sgm1_mc_1)

    L12_AVG_E1_MC_2.append(Lum_lne12_avg_e_l_sgm1_mc_2)
    L12_AVG_E2_MC_2.append(Lum_lne12_avg_e_s_sgm1_mc_2)
    L12_MED_E1_MC_2.append(Lum_lne12_med_e_l_sgm1_mc_2)
    L12_MED_E2_MC_2.append(Lum_lne12_med_e_s_sgm1_mc_2)
    PLT_ARR = [
                Z1_AVG_E1_MC,Z1_AVG_E2_MC, #1
                Z2_AVG_E1_MC,Z2_AVG_E2_MC, #3
                Z1_MED_E1_MC,Z1_MED_E2_MC, #5
                Z2_MED_E1_MC,Z2_MED_E2_MC, #7
                V1_AVG_E1_MC,V1_AVG_E2_MC, #9
                V2_AVG_E1_MC,V2_AVG_E2_MC, #11
                V1_MED_E1_MC,V1_MED_E2_MC, #13
                V2_MED_E1_MC,V2_MED_E2_MC, #15
                S1_AVG_E1_MC,S1_AVG_E2_MC, #17
                S2_AVG_E1_MC,S2_AVG_E2_MC, #19
                S1_MED_E1_MC,S1_MED_E2_MC, #21
                S2_MED_E1_MC,S2_MED_E2_MC, #23
                L1_AVG_E1_MC_1,L1_AVG_E2_MC_1, #25
                L1_MED_E1_MC_1,L1_MED_E2_MC_1, #27
                L1_AVG_E1_MC_2,L1_AVG_E2_MC_2, #29
                L1_MED_E1_MC_2,L1_MED_E2_MC_2, #31
                L2_AVG_E1_MC_1,L2_AVG_E2_MC_1, #33
                L2_MED_E1_MC_1,L2_MED_E2_MC_1, #35
                L2_AVG_E1_MC_2,L2_AVG_E2_MC_2, #37
                L2_MED_E1_MC_2,L2_MED_E2_MC_2, #39
                S12_AVG_E1_MC,S12_AVG_E2_MC,   #41
                S12_MED_E1_MC,S12_MED_E2_MC,   #43
                L12_AVG_E1_MC_1,L12_AVG_E2_MC_1, #45
                L12_MED_E1_MC_1,L12_MED_E2_MC_1, #47
                L12_AVG_E1_MC_2,L12_AVG_E2_MC_2, #49
                L12_MED_E1_MC_2,L12_MED_E2_MC_2  #51
            ]
    ################PLOT-ARRAYS######################

    ###################TABLES########################
    SMPL.append(str(nmb_smpls_mcmc))

    ZA1_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,rds_hdr))
    ZA1_MSR_MC_1_SG_L.append(z_lne1_avg_e_l_sgm1_mc)
    ZA1_MSR_MC_1_SG_H.append(z_lne1_avg_e_s_sgm1_mc)
    ZA1_MSR_MC_2_SG_L.append(z_lne1_avg_e_l_sgm2_mc)
    ZA1_MSR_MC_2_SG_H.append(z_lne1_avg_e_s_sgm2_mc)
    ZA1_MSR_MC_3_SG_L.append(z_lne1_avg_e_l_sgm3_mc)
    ZA1_MSR_MC_3_SG_H.append(z_lne1_avg_e_s_sgm3_mc)

    ZA2_MSR.append(Header_Get(spec_file_plt_lne2_avg_mcmc,rds_hdr))
    ZA2_MSR_MC_1_SG_L.append(z_lne2_avg_e_l_sgm1_mc)
    ZA2_MSR_MC_1_SG_H.append(z_lne2_avg_e_s_sgm1_mc)
    ZA2_MSR_MC_2_SG_L.append(z_lne2_avg_e_l_sgm2_mc)
    ZA2_MSR_MC_2_SG_H.append(z_lne2_avg_e_s_sgm2_mc)
    ZA2_MSR_MC_3_SG_L.append(z_lne2_avg_e_l_sgm3_mc)
    ZA2_MSR_MC_3_SG_H.append(z_lne2_avg_e_s_sgm3_mc)

    ZA12_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,rds_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,rds_hdr))
    ZA12_MSR_MC_1_SG_L.append(z_lne12_avg_e_l_sgm1_mc)
    ZA12_MSR_MC_1_SG_H.append(z_lne12_avg_e_s_sgm1_mc)
    ZA12_MSR_MC_2_SG_L.append(z_lne12_avg_e_l_sgm2_mc)
    ZA12_MSR_MC_2_SG_H.append(z_lne12_avg_e_s_sgm2_mc)
    ZA12_MSR_MC_3_SG_L.append(z_lne12_avg_e_l_sgm3_mc)
    ZA12_MSR_MC_3_SG_H.append(z_lne12_avg_e_s_sgm3_mc)

    ZM1_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,rds_hdr))
    ZM1_MSR_MC_1_SG_L.append(z_lne1_med_e_l_sgm1_mc)
    ZM1_MSR_MC_1_SG_H.append(z_lne1_med_e_s_sgm1_mc)
    ZM1_MSR_MC_2_SG_L.append(z_lne1_med_e_l_sgm2_mc)
    ZM1_MSR_MC_2_SG_H.append(z_lne1_med_e_s_sgm2_mc)
    ZM1_MSR_MC_3_SG_L.append(z_lne1_med_e_l_sgm3_mc)
    ZM1_MSR_MC_3_SG_H.append(z_lne1_med_e_s_sgm3_mc)

    ZM2_MSR.append(Header_Get(spec_file_plt_lne2_avg_mcmc,rds_hdr))
    ZM2_MSR_MC_1_SG_L.append(z_lne2_med_e_l_sgm1_mc)
    ZM2_MSR_MC_1_SG_H.append(z_lne2_med_e_s_sgm1_mc)
    ZM2_MSR_MC_2_SG_L.append(z_lne2_med_e_l_sgm2_mc)
    ZM2_MSR_MC_2_SG_H.append(z_lne2_med_e_s_sgm2_mc)
    ZM2_MSR_MC_3_SG_L.append(z_lne2_med_e_l_sgm3_mc)
    ZM2_MSR_MC_3_SG_H.append(z_lne2_med_e_s_sgm3_mc)

    ZM12_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,rds_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,rds_hdr))
    ZM12_MSR_MC_1_SG_L.append(z_lne12_med_e_l_sgm1_mc)
    ZM12_MSR_MC_1_SG_H.append(z_lne12_med_e_s_sgm1_mc)
    ZM12_MSR_MC_2_SG_L.append(z_lne12_med_e_l_sgm2_mc)
    ZM12_MSR_MC_2_SG_H.append(z_lne12_med_e_s_sgm2_mc)
    ZM12_MSR_MC_3_SG_L.append(z_lne12_med_e_l_sgm3_mc)
    ZM12_MSR_MC_3_SG_H.append(z_lne12_med_e_s_sgm3_mc)

    VA1_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,vrx_hdr))
    VA1_MSR_MC_1_SG_L.append(v_lne1_avg_e_l_sgm1_mc)
    VA1_MSR_MC_1_SG_H.append(v_lne1_avg_e_s_sgm1_mc)
    VA1_MSR_MC_2_SG_L.append(v_lne1_avg_e_l_sgm2_mc)
    VA1_MSR_MC_2_SG_H.append(v_lne1_avg_e_s_sgm2_mc)
    VA1_MSR_MC_3_SG_L.append(v_lne1_avg_e_l_sgm3_mc)
    VA1_MSR_MC_3_SG_H.append(v_lne1_avg_e_s_sgm3_mc)

    VA2_MSR.append(Header_Get(spec_file_plt_lne2_avg_mcmc,vrx_hdr))
    VA2_MSR_MC_1_SG_L.append(v_lne2_avg_e_l_sgm1_mc)
    VA2_MSR_MC_1_SG_H.append(v_lne2_avg_e_s_sgm1_mc)
    VA2_MSR_MC_2_SG_L.append(v_lne2_avg_e_l_sgm2_mc)
    VA2_MSR_MC_2_SG_H.append(v_lne2_avg_e_s_sgm2_mc)
    VA2_MSR_MC_3_SG_L.append(v_lne2_avg_e_l_sgm3_mc)
    VA2_MSR_MC_3_SG_H.append(v_lne2_avg_e_s_sgm3_mc)

    VA12_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,vrx_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,vrx_hdr))
    VA12_MSR_MC_1_SG_L.append(v_lne12_avg_e_l_sgm1_mc)
    VA12_MSR_MC_1_SG_H.append(v_lne12_avg_e_s_sgm1_mc)
    VA12_MSR_MC_2_SG_L.append(v_lne12_avg_e_l_sgm2_mc)
    VA12_MSR_MC_2_SG_H.append(v_lne12_avg_e_s_sgm2_mc)
    VA12_MSR_MC_3_SG_L.append(v_lne12_avg_e_l_sgm3_mc)
    VA12_MSR_MC_3_SG_H.append(v_lne12_avg_e_s_sgm3_mc)

    VM1_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,vrx_hdr))
    VM1_MSR_MC_1_SG_L.append(v_lne1_med_e_l_sgm1_mc)
    VM1_MSR_MC_1_SG_H.append(v_lne1_med_e_s_sgm1_mc)
    VM1_MSR_MC_2_SG_L.append(v_lne1_med_e_l_sgm2_mc)
    VM1_MSR_MC_2_SG_H.append(v_lne1_med_e_s_sgm2_mc)
    VM1_MSR_MC_3_SG_L.append(v_lne1_med_e_l_sgm3_mc)
    VM1_MSR_MC_3_SG_H.append(v_lne1_med_e_s_sgm3_mc)

    VM2_MSR.append(Header_Get(spec_file_plt_lne2_avg_mcmc,vrx_hdr))
    VM2_MSR_MC_1_SG_L.append(v_lne2_med_e_l_sgm1_mc)
    VM2_MSR_MC_1_SG_H.append(v_lne2_med_e_s_sgm1_mc)
    VM2_MSR_MC_2_SG_L.append(v_lne2_med_e_l_sgm2_mc)
    VM2_MSR_MC_2_SG_H.append(v_lne2_med_e_s_sgm2_mc)
    VM2_MSR_MC_3_SG_L.append(v_lne2_med_e_l_sgm3_mc)
    VM2_MSR_MC_3_SG_H.append(v_lne2_med_e_s_sgm3_mc)

    VM12_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,vrx_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,vrx_hdr))
    VM12_MSR_MC_1_SG_L.append(v_lne12_med_e_l_sgm1_mc)
    VM12_MSR_MC_1_SG_H.append(v_lne12_med_e_s_sgm1_mc)
    VM12_MSR_MC_2_SG_L.append(v_lne12_med_e_l_sgm2_mc)
    VM12_MSR_MC_2_SG_H.append(v_lne12_med_e_s_sgm2_mc)
    VM12_MSR_MC_3_SG_L.append(v_lne12_med_e_l_sgm3_mc)
    VM12_MSR_MC_3_SG_H.append(v_lne12_med_e_s_sgm3_mc)

    FA1_MSR.append(mu_lne1_avg)
    FA1_MSR_MC_1_SG_L.append(s_lne1_avg_e_l_sgm1_mc)
    FA1_MSR_MC_1_SG_H.append(s_lne1_avg_e_s_sgm1_mc)
    FA1_MSR_MC_2_SG_L.append(s_lne1_avg_e_l_sgm2_mc)
    FA1_MSR_MC_2_SG_H.append(s_lne1_avg_e_s_sgm2_mc)
    FA1_MSR_MC_3_SG_L.append(s_lne1_avg_e_l_sgm3_mc)
    FA1_MSR_MC_3_SG_H.append(s_lne1_avg_e_s_sgm3_mc)

    FA2_MSR.append(mu_lne2_avg)
    FA2_MSR_MC_1_SG_L.append(s_lne2_avg_e_l_sgm1_mc)
    FA2_MSR_MC_1_SG_H.append(s_lne2_avg_e_s_sgm1_mc)
    FA2_MSR_MC_2_SG_L.append(s_lne2_avg_e_l_sgm2_mc)
    FA2_MSR_MC_2_SG_H.append(s_lne2_avg_e_s_sgm2_mc)
    FA2_MSR_MC_3_SG_L.append(s_lne2_avg_e_l_sgm3_mc)
    FA2_MSR_MC_3_SG_H.append(s_lne2_avg_e_s_sgm3_mc)

    FA12_MSR.append(mu_lne1_avg/mu_lne2_avg)
    FA12_MSR_MC_1_SG_L.append(s_lne12_avg_e_l_sgm1_mc)
    FA12_MSR_MC_1_SG_H.append(s_lne12_avg_e_s_sgm1_mc)
    FA12_MSR_MC_2_SG_L.append(s_lne12_avg_e_l_sgm2_mc)
    FA12_MSR_MC_2_SG_H.append(s_lne12_avg_e_s_sgm2_mc)
    FA12_MSR_MC_3_SG_L.append(s_lne12_avg_e_l_sgm3_mc)
    FA12_MSR_MC_3_SG_H.append(s_lne12_avg_e_s_sgm3_mc)

    FM1_MSR.append(mu_lne1_med)
    FM1_MSR_MC_1_SG_L.append(s_lne1_med_e_l_sgm1_mc)
    FM1_MSR_MC_1_SG_H.append(s_lne1_med_e_s_sgm1_mc)
    FM1_MSR_MC_2_SG_L.append(s_lne1_med_e_l_sgm2_mc)
    FM1_MSR_MC_2_SG_H.append(s_lne1_med_e_s_sgm2_mc)
    FM1_MSR_MC_3_SG_L.append(s_lne1_med_e_l_sgm3_mc)
    FM1_MSR_MC_3_SG_H.append(s_lne1_med_e_s_sgm3_mc)

    FM2_MSR.append(mu_lne2_med)
    FM2_MSR_MC_1_SG_L.append(s_lne2_med_e_l_sgm1_mc)
    FM2_MSR_MC_1_SG_H.append(s_lne2_med_e_s_sgm1_mc)
    FM2_MSR_MC_2_SG_L.append(s_lne2_med_e_l_sgm2_mc)
    FM2_MSR_MC_2_SG_H.append(s_lne2_med_e_s_sgm2_mc)
    FM2_MSR_MC_3_SG_L.append(s_lne2_med_e_l_sgm3_mc)
    FM2_MSR_MC_3_SG_H.append(s_lne2_med_e_s_sgm3_mc)

    FM12_MSR.append(mu_lne1_med/mu_lne2_med)
    FM12_MSR_MC_1_SG_L.append(s_lne12_med_e_l_sgm1_mc)
    FM12_MSR_MC_1_SG_H.append(s_lne12_med_e_s_sgm1_mc)
    FM12_MSR_MC_2_SG_L.append(s_lne12_med_e_l_sgm2_mc)
    FM12_MSR_MC_2_SG_H.append(s_lne12_med_e_s_sgm2_mc)
    FM12_MSR_MC_3_SG_L.append(s_lne12_med_e_l_sgm3_mc)
    FM12_MSR_MC_3_SG_H.append(s_lne12_med_e_s_sgm3_mc)

    LA1_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,lum_hdr))
    LA1_MSR_MC_1_SG_L.append(Lum_lne1_avg_e_l_sgm1_mc_1)
    LA1_MSR_MC_1_SG_H.append(Lum_lne1_avg_e_s_sgm1_mc_1)
    LA1_MSR_MC_2_SG_L.append(Lum_lne1_avg_e_l_sgm2_mc_1)
    LA1_MSR_MC_2_SG_H.append(Lum_lne1_avg_e_s_sgm2_mc_1)
    LA1_MSR_MC_3_SG_L.append(Lum_lne1_avg_e_l_sgm3_mc_1)
    LA1_MSR_MC_3_SG_H.append(Lum_lne1_avg_e_s_sgm3_mc_1)

    LA2_MSR.append(Header_Get(spec_file_plt_lne2_avg_mcmc,lum_hdr))
    LA2_MSR_MC_1_SG_L.append(Lum_lne2_avg_e_l_sgm1_mc_1)
    LA2_MSR_MC_1_SG_H.append(Lum_lne2_avg_e_s_sgm1_mc_1)
    LA2_MSR_MC_2_SG_L.append(Lum_lne2_avg_e_l_sgm2_mc_1)
    LA2_MSR_MC_2_SG_H.append(Lum_lne2_avg_e_s_sgm2_mc_1)
    LA2_MSR_MC_3_SG_L.append(Lum_lne2_avg_e_l_sgm3_mc_1)
    LA2_MSR_MC_3_SG_H.append(Lum_lne2_avg_e_s_sgm3_mc_1)

    LA12_MSR.append(Header_Get(spec_file_plt_lne1_avg_mcmc,lum_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,lum_hdr))
    LA12_MSR_MC_1_SG_L.append(Lum_lne12_avg_e_l_sgm1_mc_1)
    LA12_MSR_MC_1_SG_H.append(Lum_lne12_avg_e_s_sgm1_mc_1)
    LA12_MSR_MC_2_SG_L.append(Lum_lne12_avg_e_l_sgm2_mc_1)
    LA12_MSR_MC_2_SG_H.append(Lum_lne12_avg_e_s_sgm2_mc_1)
    LA12_MSR_MC_3_SG_L.append(Lum_lne12_avg_e_l_sgm3_mc_1)
    LA12_MSR_MC_3_SG_H.append(Lum_lne12_avg_e_s_sgm3_mc_1)

    LM1_MSR.append(Header_Get(spec_file_plt_lne1_med_mcmc,lum_hdr))
    LM1_MSR_MC_1_SG_L.append(Lum_lne1_med_e_l_sgm1_mc_1)
    LM1_MSR_MC_1_SG_H.append(Lum_lne1_med_e_s_sgm1_mc_1)
    LM1_MSR_MC_2_SG_L.append(Lum_lne1_med_e_l_sgm2_mc_1)
    LM1_MSR_MC_2_SG_H.append(Lum_lne1_med_e_s_sgm2_mc_1)
    LM1_MSR_MC_3_SG_L.append(Lum_lne1_med_e_l_sgm3_mc_1)
    LM1_MSR_MC_3_SG_H.append(Lum_lne1_med_e_s_sgm3_mc_1)

    LM2_MSR.append(Header_Get(spec_file_plt_lne2_med_mcmc,lum_hdr))
    LM2_MSR_MC_1_SG_L.append(Lum_lne2_med_e_l_sgm1_mc_1)
    LM2_MSR_MC_1_SG_H.append(Lum_lne2_med_e_s_sgm1_mc_1)
    LM2_MSR_MC_2_SG_L.append(Lum_lne2_med_e_l_sgm2_mc_1)
    LM2_MSR_MC_2_SG_H.append(Lum_lne2_med_e_s_sgm2_mc_1)
    LM2_MSR_MC_3_SG_L.append(Lum_lne2_med_e_l_sgm3_mc_1)
    LM2_MSR_MC_3_SG_H.append(Lum_lne2_med_e_s_sgm3_mc_1)

    LM12_MSR.append(Header_Get(spec_file_plt_lne1_med_mcmc,lum_hdr)/Header_Get(spec_file_plt_lne2_med_mcmc,lum_hdr))
    LM12_MSR_MC_1_SG_L.append(Lum_lne12_med_e_l_sgm1_mc_1)
    LM12_MSR_MC_1_SG_H.append(Lum_lne12_med_e_s_sgm1_mc_1)
    LM12_MSR_MC_2_SG_L.append(Lum_lne12_med_e_l_sgm2_mc_1)
    LM12_MSR_MC_2_SG_H.append(Lum_lne12_med_e_s_sgm2_mc_1)
    LM12_MSR_MC_3_SG_L.append(Lum_lne12_med_e_l_sgm3_mc_1)
    LM12_MSR_MC_3_SG_H.append(Lum_lne12_med_e_s_sgm3_mc_1)

    LLA1_MSR.append(np.log10(Header_Get(spec_file_plt_lne1_avg_mcmc,lum_hdr)))
    LLA1_MSR_MC_1_SG_L.append(Lum_lne1_avg_e_l_sgm1_mc_2)
    LLA1_MSR_MC_1_SG_H.append(Lum_lne1_avg_e_s_sgm1_mc_2)
    LLA1_MSR_MC_2_SG_L.append(Lum_lne1_avg_e_l_sgm2_mc_2)
    LLA1_MSR_MC_2_SG_H.append(Lum_lne1_avg_e_s_sgm2_mc_2)
    LLA1_MSR_MC_3_SG_L.append(Lum_lne1_avg_e_l_sgm3_mc_2)
    LLA1_MSR_MC_3_SG_H.append(Lum_lne1_avg_e_s_sgm3_mc_2)

    LLA2_MSR.append(np.log10(Header_Get(spec_file_plt_lne2_avg_mcmc,lum_hdr)))
    LLA2_MSR_MC_1_SG_L.append(Lum_lne2_avg_e_l_sgm1_mc_2)
    LLA2_MSR_MC_1_SG_H.append(Lum_lne2_avg_e_s_sgm1_mc_2)
    LLA2_MSR_MC_2_SG_L.append(Lum_lne2_avg_e_l_sgm2_mc_2)
    LLA2_MSR_MC_2_SG_H.append(Lum_lne2_avg_e_s_sgm2_mc_2)
    LLA2_MSR_MC_3_SG_L.append(Lum_lne2_avg_e_l_sgm3_mc_2)
    LLA2_MSR_MC_3_SG_H.append(Lum_lne2_avg_e_s_sgm3_mc_2)

    LLA12_MSR.append(np.log10(Header_Get(spec_file_plt_lne1_avg_mcmc,lum_hdr)/Header_Get(spec_file_plt_lne2_avg_mcmc,lum_hdr)))
    LLA12_MSR_MC_1_SG_L.append(Lum_lne12_avg_e_l_sgm1_mc_2)
    LLA12_MSR_MC_1_SG_H.append(Lum_lne12_avg_e_s_sgm1_mc_2)
    LLA12_MSR_MC_2_SG_L.append(Lum_lne12_avg_e_l_sgm2_mc_2)
    LLA12_MSR_MC_2_SG_H.append(Lum_lne12_avg_e_s_sgm2_mc_2)
    LLA12_MSR_MC_3_SG_L.append(Lum_lne12_avg_e_l_sgm3_mc_2)
    LLA12_MSR_MC_3_SG_H.append(Lum_lne12_avg_e_s_sgm3_mc_2)

    LLM1_MSR.append(np.log10(Header_Get(spec_file_plt_lne1_med_mcmc,lum_hdr)))
    LLM1_MSR_MC_1_SG_L.append(Lum_lne1_med_e_l_sgm1_mc_2)
    LLM1_MSR_MC_1_SG_H.append(Lum_lne1_med_e_s_sgm1_mc_2)
    LLM1_MSR_MC_2_SG_L.append(Lum_lne1_med_e_l_sgm2_mc_2)
    LLM1_MSR_MC_2_SG_H.append(Lum_lne1_med_e_s_sgm2_mc_2)
    LLM1_MSR_MC_3_SG_L.append(Lum_lne1_med_e_l_sgm3_mc_2)
    LLM1_MSR_MC_3_SG_H.append(Lum_lne1_med_e_s_sgm3_mc_2)

    LLM2_MSR.append(np.log10(Header_Get(spec_file_plt_lne2_med_mcmc,lum_hdr)))
    LLM2_MSR_MC_1_SG_L.append(Lum_lne2_med_e_l_sgm1_mc_2)
    LLM2_MSR_MC_1_SG_H.append(Lum_lne2_med_e_s_sgm1_mc_2)
    LLM2_MSR_MC_2_SG_L.append(Lum_lne2_med_e_l_sgm2_mc_2)
    LLM2_MSR_MC_2_SG_H.append(Lum_lne2_med_e_s_sgm2_mc_2)
    LLM2_MSR_MC_3_SG_L.append(Lum_lne2_med_e_l_sgm3_mc_2)
    LLM2_MSR_MC_3_SG_H.append(Lum_lne2_med_e_s_sgm3_mc_2)

    LLM12_MSR.append(np.log10(Header_Get(spec_file_plt_lne1_med_mcmc,lum_hdr)/Header_Get(spec_file_plt_lne2_med_mcmc,lum_hdr)))
    LLM12_MSR_MC_1_SG_L.append(Lum_lne12_med_e_l_sgm1_mc_2)
    LLM12_MSR_MC_1_SG_H.append(Lum_lne12_med_e_s_sgm1_mc_2)
    LLM12_MSR_MC_2_SG_L.append(Lum_lne12_med_e_l_sgm2_mc_2)
    LLM12_MSR_MC_2_SG_H.append(Lum_lne12_med_e_s_sgm2_mc_2)
    LLM12_MSR_MC_3_SG_L.append(Lum_lne12_med_e_l_sgm3_mc_2)
    LLM12_MSR_MC_3_SG_H.append(Lum_lne12_med_e_s_sgm3_mc_2)
    ###################TABLES########################

    #####################TABLE-DISTRIBUTIONS-LINE1-LINE2#####################
    mct01_dist = aptbl.Table()
    mct01_dist['z_z_'+label1+'_mc'] = z_z_lne1_mc
    mct01_dist['z_z_'+label2+'_mc'] = z_z_lne2_mc
    mct01_dist['v_z_'+label1+'_mc'] = v_z_lne1_mc
    mct01_dist['v_z_'+label2+'_mc'] = v_z_lne2_mc
    mct01_dist['s_'+label1+'_avg_mc'] = s_lne1_avg_mc
    mct01_dist['s_'+label1+'_med_mc'] = s_lne1_med_mc
    mct01_dist['s_'+label2+'_avg_mc'] = s_lne2_avg_mc
    mct01_dist['s_'+label2+'_med_mc'] = s_lne2_med_mc
    mct01_dist['Lum_'+label1+'_avg_mc_1'] = Lum_lne1_avg_mc_1
    mct01_dist['Lum_'+label1+'_med_mc_1'] = Lum_lne1_med_mc_1
    mct01_dist['Lum_'+label2+'_avg_mc_1'] = Lum_lne2_avg_mc_1
    mct01_dist['Lum_'+label2+'_med_mc_1'] = Lum_lne2_med_mc_1   
    mct01_dist['Lum_'+label1+'_avg_mc_2'] = Lum_lne1_avg_mc_2
    mct01_dist['Lum_'+label1+'_med_mc_2'] = Lum_lne1_med_mc_2
    mct01_dist['Lum_'+label2+'_avg_mc_2'] = Lum_lne2_avg_mc_2
    mct01_dist['Lum_'+label2+'_med_mc_2'] = Lum_lne2_med_mc_2
    mct01_dist['z_lne12_avg_mc'] = z_lne12_avg_mc
    mct01_dist['z_lne12_med_mc'] = z_lne12_med_mc
    mct01_dist['v_lne12_avg_mc'] = v_lne12_avg_mc
    mct01_dist['v_lne12_med_mc'] = v_lne12_med_mc
    mct01_dist['Lum_lne12_avg_mc_1'] = Lum_lne12_avg_mc_1
    mct01_dist['Lum_lne12_med_mc_1'] = Lum_lne12_med_mc_1

    TABLESTATNAME_1_1  = tbl_dir_res + 'Sources'+ '-MC-'+str(iterations_mc)+'-''HATLAS-12CO-13CO-'+str(var_smpls_mcmc) +'-'+str(nmb_smpls_mcmc)+'-M'+str(method)+'.dat'
    TABLESTATNAME_1_2  = tbl_dir_res + 'Sources'+ '-MC-'+str(iterations_mc)+'-''HATLAS-12CO-13CO-'+str(var_smpls_mcmc) +'-'+str(nmb_smpls_mcmc)+'-M'+str(method)+'.csv'
    mct01_dist.write(TABLESTATNAME_1_1, format='ascii.fixed_width_two_line', overwrite = True)
    mct01_dist.write(TABLESTATNAME_1_2, format=tbl_format_opt, overwrite = True)

    print
    print (colored('Table with generated distributions (z, flux, luminosity): '+TABLESTATNAME_1_1,'green'))
    print (colored('Table with generated distributions (z, flux, luminosity): '+TABLESTATNAME_1_2,'green'))
    print
    #####################TABLE-DISTRIBUTIONS-LINE1-LINE2#####################

    print ('Sources'+ '-MC-'+str(iterations_mc)+'-''HATLAS-12CO-13CO-'+str(sbsms) +'-'+str(nmb_smpls_mcmc)+'-M'+str(method)+'.csv'  )
    if label1 == '12CO':
        label1_plt = '$^{12}\mathrm{CO}$'
    elif label1 == '13CO':
        label1_plt = '$^{13}\mathrm{CO}$'
    elif label1 == '18CO':
        label1_plt = '$\mathrm{C}^{18}\mathrm{O}$'
    elif label1 == 'LFIR':
        label1_plt = '$\mathrm{L}_{\mathrm{FIR}}$'
    else:
        pass

    if label2 == '12CO':
        label2_plt = '$^{12}\mathrm{CO}$'
    elif label2 == '13CO':
        label2_plt = '$^{13}\mathrm{CO}$'
    elif label2 == '18CO':
        label2_plt = '$\mathrm{C}^{18}\mathrm{O}$'
    elif label2 == 'LFIR':
        label2_plt = '$\mathrm{L}_{\mathrm{FIR}}$'
    else:
        pass        
 
    if label3 == '12CO':
        label3_plt = '$^{12}\mathrm{CO}$'
    elif label3 == '13CO':
        label3_plt = '$^{13}\mathrm{CO}$'
    elif label3 == '18CO':
        label3_plt = '$\mathrm{C}^{18}\mathrm{O}$'
    elif label3 == 'LFIR':
        label3_plt = '$\mathrm{L}_{\mathrm{FIR}}$'
    else:
        pass

    if plot_dist_hist == True:
        #################################################Plot##################################################
        #var_smpls_mcmc nmb_smpls_mcmc
        if dest_dir_plt != None:
            PLOTFILENAME  = str(dest_dir_plt) +  'Sources'+ '-MC-'+str(iterations_mc)+'-''HATLAS-'+label1+'-' + label2+'-'+str(var_smpls_mcmc) +'-'+str(nmb_smpls_mcmc)+'-M'+str(method)+'.pdf'
        elif dest_dir_plt == None:
            PLOTFILENAME  = mcm_dir_plt       +  'Sources'+ '-MC-'+str(iterations_mc)+'-''HATLAS-'+label1+'-' + label2+'-'+str(var_smpls_mcmc) +'-'+str(nmb_smpls_mcmc)+'-M'+str(method)+'.pdf'
        fig, ax = plt.subplots(figsize=(8, 4))

        fxsize=11*2
        fysize=8*6
        f = plt.figure(num=None, figsize=(11, 8), dpi=180, facecolor='w',
        edgecolor='k')
        plt.subplots_adjust(
        left    = (46/25.4)/fxsize,         #22-def --> 26 bigger
        bottom  = (100/25.4)/fysize,        #19-def --> 20 bigger
        right   = 1 - (18/25.4)/fxsize,     # 2-def --> 6  bigger 
        top     = 1 - (18/25.4)/fysize)     # 4-def --> 8  bigger
        plt.subplots_adjust(hspace=0.8,wspace=0.3)

        gs0 = gridspec.GridSpec(1, 1)

        ################################################Plot 1#################################################
        ########################################PREVIOUS DISTRIBUTION 1########################################
        gs11 = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs0[0])
            
        ax110 = plt.Subplot(f, gs11[0,0])
        f.add_subplot(ax110)

        ax110.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax110.xaxis.set_tick_params(labelsize=20)
        ax110.yaxis.set_tick_params(labelsize=20)
        xticklabels = ax110.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax110.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax110.minorticks_on()

        plt.xlabel('$\lambda$',fontsize=16)
        plt.ylabel('F$_\lambda$ (Jy)',fontsize=16)
        ax110.set_xlabel('z')
        ax110.set_ylabel('N/N$_{T}$')

        ax110.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax110.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)


        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax110.hist(z_1, n_bins, density=1, histtype='step',cumulative=True,#stepfilled
                                        color='blue',alpha=0.7,label=label1_plt + ' Observed')
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax110.hist(z_2, n_bins, density=1, histtype='step',cumulative=True,#stepfilled
                                        color='red',alpha=0.7,label=label2_plt + ' Observed')

        ax110.grid(True)

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)
        ########################################PREVIOUS DISTRIBUTION 1########################################
        ################################################Plot 1#################################################


        ################################################Plot 2#################################################
        ########################################PREVIOUS DISTRIBUTION 2########################################     
        ax120 = plt.Subplot(f, gs11[0,1])
        f.add_subplot(ax120)

        ax120.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax120.xaxis.set_tick_params(labelsize=16)
        ax120.yaxis.set_tick_params(labelsize=16)
        xticklabels = ax120.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax120.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax120.minorticks_on()

        plt.xlabel('$\lambda$',fontsize=16)
        plt.ylabel('F$_\lambda$ (Jy)',fontsize=16)
        ax120.set_xlabel(Splt_Hdr_Cmt)
        ax120.set_ylabel('N/N$_{T}$')

        ax120.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax120.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        # plot the cumulative histogram
        n_v_lne1, bins_v_lne1, patches_v_lne1 = ax120.hist(v_1, n_bins, density=1, histtype='step',cumulative=True,
                                        color='blue',alpha=0.7,label=label1_plt + ' Observed')
        n_v_lne2, bins_v_lne2, patches_v_lne2 = ax120.hist(v_2, n_bins, density=1, histtype='step',cumulative=True,
                                        color='red',alpha=0.7,label=label2_plt + ' Observed')

        ax120.grid(True)

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)
        ########################################PREVIOUS DISTRIBUTION 2########################################     
        ################################################Plot 2#################################################

        n_bins = int(iterations_mc / 40 )#90

        ################################################Plot 3#################################################
        #############################################MC DISTRIBUTION 1#########################################
        ax130 = plt.Subplot(f, gs11[1,0])
        f.add_subplot(ax130)

        ax130.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax130.xaxis.set_tick_params(labelsize=16)
        ax130.yaxis.set_tick_params(labelsize=16)
        xticklabels = ax130.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax130.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax130.minorticks_on()

        plt.xlabel('z',fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax130.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax130.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax130.hist(z_z_lne1_mc, n_bins, density=False, histtype='step',cumulative=True,#stepfilled
                                        color='blue',alpha=0.7,label=label1_plt + ' MCMC')
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax130.hist(z_z_lne2_mc, n_bins, density=False, histtype='step',cumulative=True,#stepfilled
                                        color='red',alpha=0.7,label=label2_plt + ' MCMC')

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)

        #############################################MC DISTRIBUTION 1#########################################
        ################################################Plot 3#################################################
        
        ################################################Plot 4#################################################
        #############################################MC DISTRIBUTION 2#########################################
        ax140 = plt.Subplot(f, gs11[1,1])
        f.add_subplot(ax140)

        ax140.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax140.xaxis.set_tick_params(labelsize=16)
        ax140.yaxis.set_tick_params(labelsize=16)
        xticklabels = ax140.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax140.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax140.minorticks_on()

        plt.xlabel(Splt_Hdr_Cmt,fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax140.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax140.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax140.hist(v_z_lne1_mc, n_bins, density=False, histtype='step',cumulative=True,#stepfilled
                                        color='blue',alpha=0.7,label=label1_plt + ' MCMC')
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax140.hist(v_z_lne2_mc, n_bins, density=False, histtype='step',cumulative=True,#stepfilled
                                        color='red',alpha=0.7,label=label2_plt + ' MCMC')

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)
        #############################################MC DISTRIBUTION 2#########################################
        ################################################Plot 4#################################################

        ################################################Plot 5#################################################
        ##########################################S MC DISTRIBUTION 1,2########################################

        ax150 = plt.Subplot(f, gs11[2,0])
        f.add_subplot(ax150)

        ax150.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax150.xaxis.set_tick_params(labelsize=16)
        ax150.yaxis.set_tick_params(labelsize=16)
        xticklabels = ax150.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax150.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax150.minorticks_on()

        plt.xlabel('S',fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax150.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax150.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax150.hist(s_lne1_avg_mc, n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='blue',alpha=0.7,label=label1_plt+' avg',linewidth=0.5)
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax150.hist(s_lne1_med_mc, n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='black',alpha=0.7,label=label1 +' med',ls=':',linewidth=0.5)
        n_v_lne1, bins_v_lne1, patches_v_lne1 = ax150.hist(s_lne2_avg_mc, n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='red',alpha=0.7,label=label2_plt+' avg',linewidth=0.5)
        n_v_lne2, bins_v_lne2, patches_v_lne2 = ax150.hist(s_lne2_med_mc, n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='black',alpha=0.7,label=label2_plt +' med',ls='-.',linewidth=0.5)

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)
        ##########################################S MC DISTRIBUTION 1,2########################################
        ################################################Plot 5#################################################

        ################################################Plot 6#################################################
        ##########################################L MC DISTRIBUTION 1,2########################################
        ax160 = plt.Subplot(f, gs11[2,1])
        f.add_subplot(ax160)

        ax160.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax160.xaxis.set_tick_params(labelsize=16)
        ax160.yaxis.set_tick_params(labelsize=16)
        #ax160.set_title('Teste')
        xticklabels = ax160.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax160.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax160.minorticks_on()

        plt.xlabel('log$_{10}$[L]',fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax160.yaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax160.xaxis.set_tick_params(which='both',labelsize=16,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        weights = np.ones_like(~np.isnan(Lum_lne1_avg_mc_2))/float(len(~np.isnan(Lum_lne1_avg_mc_2)))
        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax160.hist(Lum_lne1_avg_mc_2[~np.isnan(Lum_lne1_avg_mc_2)], n_bins, density=False, histtype='step',cumulative=False,#stepfilled weights=weights,
                                        color='blue',alpha=0.7,label=label1_plt +' avg',linewidth=0.5)
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax160.hist(Lum_lne1_med_mc_2[~np.isnan(Lum_lne1_med_mc_2)], n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='black',alpha=0.7,label=label1_plt + ' med',ls=':',linewidth=0.5)
        n_v_lne1, bins_v_lne1, patches_v_lne1 = ax160.hist(Lum_lne2_avg_mc_2[~np.isnan(Lum_lne2_avg_mc_2)], n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='red',alpha=0.7,label=label2 +' avg',linewidth=0.5)
        n_v_lne2, bins_v_lne2, patches_v_lne2 = ax160.hist(Lum_lne2_med_mc_2[~np.isnan(Lum_lne2_med_mc_2)], n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='black',alpha=0.7,label=label2_plt + ' med',ls='-.',linewidth=0.5)
        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)
        ##########################################L MC DISTRIBUTION 1,2########################################
        ################################################Plot 6#################################################
        
        print
        print (colored('MCMC iteration number: '+str(iterations_mc),'yellow'))
        print (colored('Bin number: '+str(nbins),'yellow'))
        print
        n_bins = iterations_mc*4

        ################################################Plot 7#################################################
        #######################################S RATIO MC DISTRIBUTION 1,2#####################################

        ax170 = plt.Subplot(f, gs11[3,0])
        f.add_subplot(ax170)

        ax170.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax170.xaxis.set_tick_params(labelsize=20)
        ax170.yaxis.set_tick_params(labelsize=20)
        xticklabels = ax170.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax170.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax170.minorticks_on()

        plt.xlabel('S $_{'+label1+'}$'+'/'+'S $_{'+label2+'}$',fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax170.yaxis.set_tick_params(which='both',labelsize=20,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax170.xaxis.set_tick_params(which='both',labelsize=20,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax170.hist(abs(s_lne1_avg_mc/s_lne2_avg_mc), n_bins, density=False, histtype='step',cumulative=False,
                                        color='blue',alpha=0.7,label=label1_plt+'/' +label2+' avg',linewidth=0.5)
        n_z_lne2, bins_z_lne2, patches_z_lne2 = ax170.hist(abs(s_lne1_med_mc/s_lne2_med_mc), n_bins, density=False, histtype='step',cumulative=False,
                                        color='black',alpha=0.7,label=label1_plt+'/' +label2 +' med',ls=':',linewidth=0.5)

        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)

        min_x =  -10
        max_x =  100
        plt.xlim([min_x,max_x])
        xmin, xmax = plt.xlim()
        plt.xlim((xmin,xmax))
        #######################################S RATIO MC DISTRIBUTION 1,2#####################################
        ################################################Plot 7#################################################

        ################################################Plot 8#################################################
        ################s#######################L RATIO MC DISTRIBUTION 1,2#################s####################
        ax180 = plt.Subplot(f, gs11[3,1])
        f.add_subplot(ax180)

        ax180.set_rasterization_zorder(1)
        plt.autoscale(enable=True, axis='y', tight=False)
        ax180.xaxis.set_tick_params(labelsize=20)
        ax180.yaxis.set_tick_params(labelsize=20)
        xticklabels = ax180.get_xticklabels()
        plt.setp(xticklabels, visible=True)
        yticklabels = ax180.get_yticklabels()
        plt.setp(yticklabels, visible=True)

        plt.tick_params(which='both', width=1.0)
        plt.tick_params(which='major', length=10)
        plt.tick_params(which='minor', length=5)
        ax180.minorticks_on()

        plt.xlabel('L $_{'+label1+'}$'+'/'+'L $_{'+label2+'}$',fontsize=16)
        plt.ylabel('N',fontsize=16)

        ax180.yaxis.set_tick_params(which='both',labelsize=20,direction='in',color='black',bottom=True,top=True,left=True,right=True)
        ax180.xaxis.set_tick_params(which='both',labelsize=20,direction='in',color='black',bottom=True,top=True,left=True,right=True)

        weights = np.ones_like(~np.isnan(Lum_lne1_avg_mc_2))/float(len(~np.isnan(Lum_lne1_avg_mc_2)))

        Lum_lne1_avg_mc_2[~np.isnan(Lum_lne2_avg_mc_2)]/Lum_lne2_avg_mc_2[~np.isnan(Lum_lne2_avg_mc_2)]
        Lum_lne1_med_mc_2[~np.isnan(Lum_lne2_med_mc_2)]/Lum_lne2_med_mc_2[~np.isnan(Lum_lne2_med_mc_2)]

        Lum_lne1_lne2_avg_mc_1_rto = Lum_lne1_avg_mc_1/Lum_lne2_avg_mc_1
        Lum_lne1_lne2_med_mc_1_rto = Lum_lne1_med_mc_1/Lum_lne2_med_mc_1
        n_z_lne1, bins_z_lne1, patches_z_lne1 = ax180.hist(Lum_lne1_lne2_avg_mc_1_rto[~np.isnan(Lum_lne1_lne2_avg_mc_1_rto)],
                                        n_bins, density=False, histtype='step',cumulative=False,#stepfilled weights=weights,
                                        color='blue',alpha=0.7,label=label1_plt+'/' +label2+' avg',linewidth=0.5)
        n_v_lne1, bins_v_lne1, patches_v_lne1 = ax180.hist(Lum_lne1_lne2_med_mc_1_rto[~np.isnan(Lum_lne1_lne2_med_mc_1_rto)],
                                        n_bins, density=False, histtype='step',cumulative=False,#stepfilled
                                        color='black',alpha=0.7,label=label1_plt+'/' +label2 +' med',ls=':',linewidth=0.5)
        lg=plt.legend(loc=1,prop={'size':12})
        lg.draw_frame(False)

        min_x =  -10
        max_x =  100
        plt.xlim([min_x,max_x])
        xmin, xmax = plt.xlim()
        plt.xlim((xmin,xmax))       
        #######################################S RATIO MC DISTRIBUTION 1,2#####################################
        ################################################Plot 8#################################################
        plt.savefig(PLOTFILENAME)
        print (colored('Generated Plot: '+str(PLOTFILENAME),'cyan'))
        plt.close('all')
        #################################################Plot##################################################
    elif plot_dist_hist == False:
        pass

    MCMC_OPT = [
                FA1_MSR,                                          #0       
                FA1_MSR_MC_1_SG_L,FA1_MSR_MC_1_SG_H,              #2
                FA1_MSR_MC_2_SG_L,FA1_MSR_MC_2_SG_H,              #4
                FA1_MSR_MC_3_SG_L,FA1_MSR_MC_3_SG_H,              #6
                FA2_MSR,                                          #7       
                FA2_MSR_MC_1_SG_L,FA2_MSR_MC_1_SG_H,              #9
                FA2_MSR_MC_2_SG_L,FA2_MSR_MC_2_SG_H,              #11
                FA2_MSR_MC_3_SG_L,FA2_MSR_MC_3_SG_H,              #13
                FA12_MSR,                                         #14       
                FA12_MSR_MC_1_SG_L,FA12_MSR_MC_1_SG_H,            #16
                FA12_MSR_MC_2_SG_L,FA12_MSR_MC_2_SG_H,            #18
                FA12_MSR_MC_3_SG_L,FA12_MSR_MC_3_SG_H,            #20
                FM1_MSR,                                          #21
                FM1_MSR_MC_1_SG_L,FM1_MSR_MC_1_SG_H,              #23
                FM1_MSR_MC_2_SG_L,FM1_MSR_MC_2_SG_H,              #25
                FM1_MSR_MC_3_SG_L,FM1_MSR_MC_3_SG_H,              #27
                FM2_MSR,                                          #28
                FM2_MSR_MC_1_SG_L,FM2_MSR_MC_1_SG_H,              #30
                FM2_MSR_MC_2_SG_L,FM2_MSR_MC_2_SG_H,              #32
                FM2_MSR_MC_3_SG_L,FM2_MSR_MC_3_SG_H,              #34
                FM12_MSR,                                         #35
                FM12_MSR_MC_1_SG_L,FM12_MSR_MC_1_SG_H,            #37
                FM12_MSR_MC_2_SG_L,FM12_MSR_MC_2_SG_H,            #39
                FM12_MSR_MC_3_SG_L,FM12_MSR_MC_3_SG_H,            #41
                LA1_MSR,                                          #42
                LA1_MSR_MC_1_SG_L,LA1_MSR_MC_1_SG_H,              #44
                LA1_MSR_MC_2_SG_L,LA1_MSR_MC_2_SG_H,              #46
                LA1_MSR_MC_3_SG_L,LA1_MSR_MC_3_SG_H,              #48
                LA2_MSR,                                          #49       
                LA2_MSR_MC_1_SG_L,LA2_MSR_MC_1_SG_H,              #51
                LA2_MSR_MC_2_SG_L,LA2_MSR_MC_2_SG_H,              #53
                LA2_MSR_MC_3_SG_L,LA2_MSR_MC_3_SG_H,              #55
                LA12_MSR,                                         #56
                LA12_MSR_MC_1_SG_L,LA12_MSR_MC_1_SG_H,            #58
                LA12_MSR_MC_2_SG_L,LA12_MSR_MC_2_SG_H,            #60
                LA12_MSR_MC_3_SG_L,LA12_MSR_MC_3_SG_H,            #62
                LM1_MSR,                                          #63       
                LM1_MSR_MC_1_SG_L,LM1_MSR_MC_1_SG_H,              #65
                LM1_MSR_MC_2_SG_L,LM1_MSR_MC_2_SG_H,              #67
                LM1_MSR_MC_3_SG_L,LM1_MSR_MC_3_SG_H,              #69
                LM2_MSR,                                          #70
                LM2_MSR_MC_1_SG_L,LM2_MSR_MC_1_SG_H,              #72
                LM2_MSR_MC_2_SG_L,LM2_MSR_MC_2_SG_H,              #74
                LM2_MSR_MC_3_SG_L,LM2_MSR_MC_3_SG_H,              #76
                LM12_MSR,                                         #77
                LM12_MSR_MC_1_SG_L,LM12_MSR_MC_1_SG_H,            #79
                LM12_MSR_MC_2_SG_L,LM12_MSR_MC_2_SG_H,            #81
                LM12_MSR_MC_3_SG_L,LM12_MSR_MC_3_SG_H,            #83
                LLA1_MSR,                                         #84
                LLA1_MSR_MC_1_SG_L,LLA1_MSR_MC_1_SG_H,            #86
                LLA1_MSR_MC_2_SG_L,LLA1_MSR_MC_2_SG_H,            #88
                LLA1_MSR_MC_3_SG_L,LLA1_MSR_MC_3_SG_H,            #90
                LLA2_MSR,                                         #91
                LLA2_MSR_MC_1_SG_L,LLA2_MSR_MC_1_SG_H,            #93
                LLA2_MSR_MC_2_SG_L,LLA2_MSR_MC_2_SG_H,            #95
                LLA2_MSR_MC_3_SG_L,LLA2_MSR_MC_3_SG_H,            #97
                LLA12_MSR,                                        #98
                LLA12_MSR_MC_1_SG_L,LLA12_MSR_MC_1_SG_H,          #100
                LLA12_MSR_MC_2_SG_L,LLA12_MSR_MC_2_SG_H,          #102
                LLA12_MSR_MC_3_SG_L,LLA12_MSR_MC_3_SG_H,          #104
                LLM1_MSR,                                         #105
                LLM1_MSR_MC_1_SG_L,LLM1_MSR_MC_1_SG_H,            #107
                LLM1_MSR_MC_2_SG_L,LLM1_MSR_MC_2_SG_H,            #109
                LLM1_MSR_MC_3_SG_L,LLM1_MSR_MC_3_SG_H,            #111
                LLM2_MSR,                                         #112
                LLM2_MSR_MC_1_SG_L,LLM2_MSR_MC_1_SG_H,            #114
                LLM2_MSR_MC_2_SG_L,LLM2_MSR_MC_2_SG_H,            #116
                LLM2_MSR_MC_3_SG_L,LLM2_MSR_MC_3_SG_H,            #118
                LLM12_MSR,                                        #119
                LLM12_MSR_MC_1_SG_L,LLM12_MSR_MC_1_SG_H,          #121
                LLM12_MSR_MC_2_SG_L,LLM12_MSR_MC_2_SG_H,          #123
                LLM12_MSR_MC_3_SG_L,LLM12_MSR_MC_3_SG_H,          #125
                VA1_MSR,                                          #126
                VA1_MSR_MC_1_SG_L,VA1_MSR_MC_1_SG_H,              #128
                VA1_MSR_MC_2_SG_L,VA1_MSR_MC_2_SG_H,              #130
                VA1_MSR_MC_3_SG_L,VA1_MSR_MC_3_SG_H,              #132
                VA2_MSR,                                          #133
                VA2_MSR_MC_1_SG_L,VA2_MSR_MC_1_SG_H,              #135
                VA2_MSR_MC_2_SG_L,VA2_MSR_MC_2_SG_H,              #137
                VA2_MSR_MC_3_SG_L,VA2_MSR_MC_3_SG_H,              #139
                VA12_MSR,                                         #140
                VA12_MSR_MC_1_SG_L,VA12_MSR_MC_1_SG_H,            #142
                VA12_MSR_MC_2_SG_L,VA12_MSR_MC_2_SG_H,            #144
                VA12_MSR_MC_3_SG_L,VA12_MSR_MC_3_SG_H,            #146
                VM1_MSR,                                          #147
                VM1_MSR_MC_1_SG_L,VM1_MSR_MC_1_SG_H,              #149
                VM1_MSR_MC_2_SG_L,VM1_MSR_MC_2_SG_H,              #151
                VM1_MSR_MC_3_SG_L,VM1_MSR_MC_3_SG_H,              #153
                VM2_MSR,                                          #154
                VM2_MSR_MC_1_SG_L,VM2_MSR_MC_1_SG_H,              #156
                VM2_MSR_MC_2_SG_L,VM2_MSR_MC_2_SG_H,              #158
                VM2_MSR_MC_3_SG_L,VM2_MSR_MC_3_SG_H,              #160
                VM12_MSR,                                         #161
                VM12_MSR_MC_1_SG_L,VM12_MSR_MC_1_SG_H,            #163
                VM12_MSR_MC_2_SG_L,VM12_MSR_MC_2_SG_H,            #165
                VM12_MSR_MC_3_SG_L,VM12_MSR_MC_3_SG_H,            #167
                ZA1_MSR,                                          #168
                ZA1_MSR_MC_1_SG_L,ZA1_MSR_MC_1_SG_H,              #170
                ZA1_MSR_MC_2_SG_L,ZA1_MSR_MC_2_SG_H,              #172
                ZA1_MSR_MC_3_SG_L,ZA1_MSR_MC_3_SG_H,              #174
                ZA2_MSR,                                          #176
                ZA2_MSR_MC_1_SG_L,ZA2_MSR_MC_1_SG_H,              #178
                ZA2_MSR_MC_2_SG_L,ZA2_MSR_MC_2_SG_H,              #180
                ZA2_MSR_MC_3_SG_L,ZA2_MSR_MC_3_SG_H,              #182
                ZA12_MSR,                                         #183
                ZA12_MSR_MC_1_SG_L,ZA12_MSR_MC_1_SG_H,            #185
                ZA12_MSR_MC_2_SG_L,ZA12_MSR_MC_2_SG_H,            #187
                ZA12_MSR_MC_3_SG_L,ZA12_MSR_MC_3_SG_H,            #199
                ZM1_MSR,                                          #191
                ZM1_MSR_MC_1_SG_L,ZM1_MSR_MC_1_SG_H,              #193
                ZM1_MSR_MC_2_SG_L,ZM1_MSR_MC_2_SG_H,              #195
                ZM1_MSR_MC_3_SG_L,ZM1_MSR_MC_3_SG_H,              #197
                ZM2_MSR,                                          #198
                ZM2_MSR_MC_1_SG_L,ZM2_MSR_MC_1_SG_H,              #200
                ZM2_MSR_MC_2_SG_L,ZM2_MSR_MC_2_SG_H,              #202
                ZM2_MSR_MC_3_SG_L,ZM2_MSR_MC_3_SG_H,              #204
                ZM12_MSR,                                         #205
                ZM12_MSR_MC_1_SG_L,ZM12_MSR_MC_1_SG_H,            #207
                ZM12_MSR_MC_2_SG_L,ZM12_MSR_MC_2_SG_H,            #208
                ZM12_MSR_MC_3_SG_L,ZM12_MSR_MC_3_SG_H             #209
            ]
    return MCMC_OPT,PLT_ARR

def MCMC_generator(iterations_mc,line1,line2,method,error,sbsms_mcmc,sbsmn_mcmc,*args,**kwargs):
    plt_scl         = kwargs.get('plt_scl',None)
    log_lm          = kwargs.get('log_lm','error')
    func1           = kwargs.get('func1','avg')
    func2           = kwargs.get('func2','med')
    plot_dist_hist  = kwargs.get('plot_dist_hist',True)
    spc_wdt_dir     = kwargs.get('spc_wdt_dir'    ,500)
    mask_radi_as_ms = kwargs.get('mask_radi_as_ms',15)

    Splt_Vars = split_variable_vars(sbsms_mcmc)
    Tbl_Splt_Col = Splt_Vars[0]
    Tbl_Splt_Hdr = Splt_Vars[1]
    Tbl_Splt_Hdr_Cmt = Splt_Vars[2]
    Splt_CNm = Splt_Vars[3]
    Splt_Hdr = Splt_Vars[1]
    Splt_Hdr_Plt = Splt_Vars[4]
    Splt_Plt_lbl = Splt_Vars[5]

    #DEFINES HEADERS TO READ FLX, LUM AND ERRS#
    if method == 1 and log_lm == None:
        #1
        #STT_TFL  TFlx SUM All Chns * CbeWth               #1
        #STT_LMT  TFlx SUM All CH X CbeWth 2 Lum                   #2
        #STT_LLT  TFlx SUM All CH X CbeWth 2 Lum [log]             #3
        #STT_TFE  TFlx SUM All Chns * CbeWth Err           #1
        #STT_SL1  TFlx SUM All CH 2 Lum Err 1 sgm lw lmt 15.9 pct  #2
        #STT_SH1  TFlx SUM All CH 2 Lum Err 1 sgm hg lmt 84.1 pct  #2
        #STT_LL1  TFlx SUM All CH 2 Lum [log] Err 1 sgm lw lmt 15  #3
        #STT_LH1  TFlx SUM All CH 2 Lum [log] Err 1 sgm hg lmt 84  #3
        flx_hdr     = 'STT_TFL'
        lum_hdr     = 'STT_LMT'
        flx_hdr_e   = 'STT_TFE'
        lum_hdr_e1  = 'STT_SL'+str(error)
        lum_hdr_e2  = 'STT_SH'+str(error)

        lum_hdr_e11 = 'STT_SL1'
        lum_hdr_e21 = 'STT_SH1'
        lum_hdr_e12 = 'STT_SL2'
        lum_hdr_e22 = 'STT_SH2'
        lum_hdr_e13 = 'STT_SL3'
        lum_hdr_e23 = 'STT_SH3'

    elif method == 1 and log_lm == 'error':
        flx_hdr     = 'STT_TFL'
        lum_hdr     = 'STT_LMT'
        flx_hdr_e   = 'STT_TFE'
        lum_hdr_e1  = 'STT_LL'+str(error)
        lum_hdr_e2  = 'STT_LH'+str(error)

        lum_hdr_e11 = 'STT_LL1'
        lum_hdr_e21 = 'STT_LH1'
        lum_hdr_e12 = 'STT_LL2'
        lum_hdr_e22 = 'STT_LH2'
        lum_hdr_e13 = 'STT_LL3'
        lum_hdr_e23 = 'STT_LH3'

    elif method == 1 and log_lm == 'value':
        flx_hdr     = 'STT_TFL'
        lum_hdr     = 'STT_LLT'
        flx_hdr_e   = 'STT_TFE'
        lum_hdr_e1  = 'STT_SL'+str(error)
        lum_hdr_e2  = 'STT_SH'+str(error)

        lum_hdr_e11 = 'STT_SL1'
        lum_hdr_e21 = 'STT_SH1'
        lum_hdr_e12 = 'STT_SL2'
        lum_hdr_e22 = 'STT_SH2'
        lum_hdr_e13 = 'STT_SL3'
        lum_hdr_e23 = 'STT_SH3'

    elif method == 1 and log_lm == 'both':
        flx_hdr     = 'STT_TFL'
        lum_hdr     = 'STT_LLT'
        flx_hdr_e   = 'STT_TFE'
        lum_hdr_e1  = 'STT_LL'+str(error)
        lum_hdr_e2  = 'STT_LH'+str(error)

        lum_hdr_e11 = 'STT_LL1'
        lum_hdr_e21 = 'STT_LH1'
        lum_hdr_e12 = 'STT_LL2'
        lum_hdr_e22 = 'STT_LH2'
        lum_hdr_e13 = 'STT_LL3'
        lum_hdr_e23 = 'STT_LH3'

    elif method == 2 and log_lm == None:
        #2
        #FTS_A2M  1DGF Area M SUM                                  #1
        #FTS_LUM  1DGF Ar2Lum SUM                                  #2
        #FTS_LLM  1DGF Ar2Lum [log] SUM                            #3
        #FTS_AME  1DGF Area M Err SUM                              #1
        #FMS_ML1  1DGF Ar2Lum Lum M SUM Err 1 sgm lw lmt 15.9 pct  #2
        #FMS_MH1  1DGF Ar2Lum Lum M SUM Err 1 sgm hg lmt 84.1 pct  #2
        #FMS_LL1  1DGF Ar2Lum Lum M [log] SUM Err 1 sgm lw lmt 15  #3
        #FMS_LH1  1DGF Ar2Lum Lum M [log] SUM Err 1 sgm hg lmt 84  #3
        flx_hdr     = 'FTS_A2M'
        lum_hdr     = 'FTS_LUM'
        flx_hdr_e   = 'FTS_AME'
        lum_hdr_e1  = 'FMS_ML'+str(error)
        lum_hdr_e2  = 'FMS_MH'+str(error)

        lum_hdr_e11 = 'FMS_ML1'
        lum_hdr_e21 = 'FMS_MH1'
        lum_hdr_e12 = 'FMS_ML2'
        lum_hdr_e22 = 'FMS_MH2'
        lum_hdr_e13 = 'FMS_ML3'
        lum_hdr_e23 = 'FMS_MH3'

    elif method == 2 and log_lm == 'error':
        flx_hdr     = 'FTS_A2M'
        lum_hdr     = 'FTS_LUM'
        flx_hdr_e   = 'FTS_AME'
        lum_hdr_e1  = 'FMS_LL'+str(error)
        lum_hdr_e2  = 'FMS_LH'+str(error)

        lum_hdr_e11 = 'FMS_LL1'
        lum_hdr_e21 = 'FMS_LH1'
        lum_hdr_e12 = 'FMS_LL2'
        lum_hdr_e22 = 'FMS_LH2'
        lum_hdr_e13 = 'FMS_LL3'
        lum_hdr_e23 = 'FMS_LH3'

    elif method == 2 and log_lm == 'value':
        flx_hdr     = 'FTS_A2M'
        lum_hdr     = 'FTS_LLM'
        flx_hdr_e   = 'FTS_AME'
        lum_hdr_e1  = 'FMS_ML'+str(error)
        lum_hdr_e2  = 'FMS_MH'+str(error)

        lum_hdr_e11 = 'FMS_ML1'
        lum_hdr_e21 = 'FMS_MH1'
        lum_hdr_e12 = 'FMS_ML2'
        lum_hdr_e22 = 'FMS_MH2'
        lum_hdr_e13 = 'FMS_ML3'
        lum_hdr_e23 = 'FMS_MH3'

    elif method == 2 and log_lm == 'both':
        flx_hdr     = 'FTS_A2M'
        lum_hdr     = 'FTS_LLM'
        flx_hdr_e   = 'FTS_AME'
        lum_hdr_e1  = 'FMS_LL'+str(error)
        lum_hdr_e2  = 'FMS_LH'+str(error)

        lum_hdr_e11 = 'FMS_LL1'
        lum_hdr_e21 = 'FMS_LH1'
        lum_hdr_e12 = 'FMS_LL2'
        lum_hdr_e22 = 'FMS_LH2'
        lum_hdr_e13 = 'FMS_LL3'
        lum_hdr_e23 = 'FMS_LH3'

    elif method == 3 and log_lm == None:
        #3
        #S2G_FT1  2DGF Vol X CbeWth SUM                         #1
        #S2G_LMT  2DGF Vol X CbeWth Lum SUM                     #3
        #S2G_LLT  2DGF Vol X CbeWth Lum [log] SUM               #4
        #S2G_F1E  2DGF Vol X CbeWth SUM Err                     #1
        #S2L_T1L  2DGF Vol Lum T X CbeWth SUM Err 1 sgm lw lm   #3
        #S2L_T1H  2DGF Vol Lum T X CbeWth SUM Err 1 sgm hg lm   #3
        #S2M_T1L  2DGF Vol Lum T X CbeWth [log] SUM Err 1 sgm   #4
        #S2M_T1H  2DGF Vol Lum T X CbeWth [log] SUM Err 1 sgm   #4
        flx_hdr     = 'S2G_FT1'#'S2G_FT1'
        lum_hdr     = 'S2G_LM1'#'S2G_LMT'
        flx_hdr_e   = 'S2G_F1E'#'S2G_F1E'
        lum_hdr_e1  = 'S2L_1'+str(error)+'L'#'S2L_T'+str(error)+'L'
        lum_hdr_e2  = 'S2L_1'+str(error)+'H'#'S2L_T'+str(error)+'H'

        lum_hdr_e11 = 'S2L_11L'#'S2L_T1L'
        lum_hdr_e21 = 'S2L_11H'#'S2L_T1H'
        lum_hdr_e12 = 'S2L_12L'#'S2L_T2L'
        lum_hdr_e22 = 'S2L_12H'#'S2L_T2H'
        lum_hdr_e13 = 'S2L_13L'#'S2L_T3L'
        lum_hdr_e23 = 'S2L_13H'#'S2L_T3H'

    elif method == 3 and log_lm == 'error':
        flx_hdr     = 'S2G_FT1'#'S2G_FT1'
        lum_hdr     = 'S2G_LM1'#'S2G_LMT'
        flx_hdr_e   = 'S2G_F1E'#'S2G_F1E'
        lum_hdr_e1  = 'S2M_1'+str(error)+'L'#'S2M_T'+str(error)+'L'
        lum_hdr_e2  = 'S2M_1'+str(error)+'H'#'S2M_T'+str(error)+'H'

        lum_hdr_e11 = 'S2M_11L'#'S2M_T1L'
        lum_hdr_e21 = 'S2M_11H'#'S2M_T1H'
        lum_hdr_e12 = 'S2M_12L'#'S2M_T2L'
        lum_hdr_e22 = 'S2M_12H'#'S2M_T2H'
        lum_hdr_e13 = 'S2M_13L'#'S2M_T3L'
        lum_hdr_e23 = 'S2M_13H'#'S2M_T3H'

    elif method == 3 and log_lm == 'value':
        flx_hdr     = 'S2G_FT1'#'S2G_FT1'
        lum_hdr     = 'S2G_LL1'#'S2G_LLT'
        flx_hdr_e   = 'S2G_F1E'#'S2G_F1E'
        lum_hdr_e1  = 'S2L_1'+str(error)+'L'#'S2L_T'+str(error)+'L'
        lum_hdr_e2  = 'S2L_1'+str(error)+'H'#'S2L_T'+str(error)+'H'

        lum_hdr_e11 = 'S2L_11L'#'S2L_T1L'
        lum_hdr_e21 = 'S2L_11H'#'S2L_T1H'
        lum_hdr_e12 = 'S2L_12L'#'S2L_T2L'
        lum_hdr_e22 = 'S2L_12H'#'S2L_T2H'
        lum_hdr_e13 = 'S2L_13L'#'S2L_T3L'
        lum_hdr_e23 = 'S2L_13H'#'S2L_T3H'

    elif method == 3 and log_lm == 'both':
        flx_hdr     = 'S2G_FT1'#'S2G_FT1'
        lum_hdr     = 'S2G_LL1'#'S2G_LLT'
        flx_hdr_e   = 'S2G_F1E'#'S2G_F1E'
        lum_hdr_e1  = 'S2M_1'+str(error)+'L'#'S2M_T'+str(error)+'L'
        lum_hdr_e2  = 'S2M_1'+str(error)+'H'#'S2M_T'+str(error)+'H'

        lum_hdr_e11 = 'S2M_11L'#'S2M_T1L'
        lum_hdr_e21 = 'S2M_11H'#'S2M_T1H'
        lum_hdr_e12 = 'S2M_12L'#'S2M_T2L'
        lum_hdr_e22 = 'S2M_12H'#'S2M_T2H'
        lum_hdr_e13 = 'S2M_13L'#'S2M_T3L'
        lum_hdr_e23 = 'S2M_13H'#'S2M_T3H'

    elif method == 4 and log_lm == None:
        #4
        #S2G_FT2  2DGF Vol X CbeWth 1DGF FWHM SUM               #2
        #S2G_LMF  2DGF Vol X CbeWth 1DGF FWHM Lum SUM           #5
        #S2G_LLF  2DGF Vol X CbeWth 1DGF FWHM Lum [log] SUM     #6
        #S2G_F2E  2DGF Vol X CbeWth 1DGF FWHM SUM Err           #2
        #S2L_F1L  2DGF Vol Lum T X CbeWth 1DGF FWHM SUM Err 1   #5
        #S2L_F1H  2DGF Vol Lum T X CbeWth 1DGF FWHM SUM Err 1   #5
        #S2M_F1L  2DGF Vol Lum T X CbeWth 1DGF FWHM [log] SUM   #6
        #S2M_F1H  2DGF Vol Lum T X CbeWth 1DGF FWHM [log] SUM   #6
        flx_hdr     = 'S2G_FT2'#'S2G_FT2'
        lum_hdr     = 'S2G_LM2'#'S2G_LMF'
        flx_hdr_e   = 'S2G_F2E'#'S2G_F2E'
        lum_hdr_e1  = 'S2L_2'+str(error)+'L'# 'S2L_F'+str(error)+'L'
        lum_hdr_e2  = 'S2L_2'+str(error)+'H'# 'S2L_F'+str(error)+'H'

        lum_hdr_e11 = 'S2L_21L'#'S2L_F1L'
        lum_hdr_e21 = 'S2L_21H'#'S2L_F1H'
        lum_hdr_e12 = 'S2L_22L'#'S2L_F2L'
        lum_hdr_e22 = 'S2L_22H'#'S2L_F2H'
        lum_hdr_e13 = 'S2L_23L'#'S2L_F3L'
        lum_hdr_e23 = 'S2L_23H'#'S2L_F3H'

    elif method == 4 and log_lm == 'error':
        flx_hdr     = 'S2G_FT2'#'S2G_FT2'
        lum_hdr     = 'S2G_LM2'#'S2G_LMF'
        flx_hdr_e   = 'S2G_F2E'#'S2G_F2E'
        lum_hdr_e1  = 'S2M_2'+str(error)+'L'#'S2M_F'+str(error)+'L'
        lum_hdr_e2  = 'S2M_2'+str(error)+'L'#'S2M_F'+str(error)+'H'

        lum_hdr_e11 = 'S2M_21L'#'S2M_F1L'
        lum_hdr_e21 = 'S2M_21H'#'S2M_F1H'
        lum_hdr_e12 = 'S2M_22L'#'S2M_F2L'
        lum_hdr_e22 = 'S2M_22H'#'S2M_F2H'
        lum_hdr_e13 = 'S2M_23L'#'S2M_F3L'
        lum_hdr_e23 = 'S2M_23H'#'S2M_F3H'

    elif method == 4 and log_lm == 'value':
        flx_hdr     = 'S2G_FT2'#'S2G_FT2'
        lum_hdr     = 'S2G_LL2'#'S2G_LLF'
        flx_hdr_e   = 'S2G_F2E'#'S2G_F2E'
        lum_hdr_e1  = 'S2L_2' + str(error) + 'L'# 'S2L_F'+str(error)+'L'
        lum_hdr_e2  = 'S2L_2' + str(error) + 'H'# 'S2L_F'+str(error)+'H'

        lum_hdr_e11 = 'S2L_21L'#'S2M_F1L'
        lum_hdr_e21 = 'S2L_21H'#'S2M_F1H'
        lum_hdr_e12 = 'S2L_22L'#'S2M_F2L'
        lum_hdr_e22 = 'S2L_22H'#'S2M_F2H'
        lum_hdr_e13 = 'S2L_23L'#'S2M_F3L'
        lum_hdr_e23 = 'S2L_23H'#'S2M_F3H'

    elif method == 4 and log_lm == 'both':
        flx_hdr     = 'S2G_FT2'#'S2G_FT2'
        lum_hdr     = 'S2G_LL2'#'S2G_LLF'
        flx_hdr_e   = 'S2G_F2E'#'S2G_F2E'
        lum_hdr_e1  = 'S2M_2' + str(error) + 'L'#'S2M_F'+str(error)+'L'
        lum_hdr_e2  = 'S2M_2' + str(error) + 'H'#'S2M_F'+str(error)+'H'

        lum_hdr_e11 = 'S2M_21L'#'S2M_F1L'
        lum_hdr_e21 = 'S2M_21H'#'S2M_F1H'
        lum_hdr_e12 = 'S2M_22L'#'S2M_F2L'
        lum_hdr_e22 = 'S2M_22H'#'S2M_F2H'
        lum_hdr_e13 = 'S2M_23L'#'S2M_F3L'
        lum_hdr_e23 = 'S2M_23H'#'S2M_F3H'

    elif method == 5 and log_lm == None:
        #4
        #S2G_FT2  2DGF Vol X CbeWth 1DGF FWHM SUM               #2
        #S2G_LMF  2DGF Vol X CbeWth 1DGF FWHM Lum SUM           #5
        #S2G_LLF  2DGF Vol X CbeWth 1DGF FWHM Lum [log] SUM     #6
        #S2G_F2E  2DGF Vol X CbeWth 1DGF FWHM SUM Err           #2
        #S2L_F1L  2DGF Vol Lum T X CbeWth 1DGF FWHM SUM Err 1   #5
        #S2L_F1H  2DGF Vol Lum T X CbeWth 1DGF FWHM SUM Err 1   #5
        #S2M_F1L  2DGF Vol Lum T X CbeWth 1DGF FWHM [log] SUM   #6
        #S2M_F1H  2DGF Vol Lum T X CbeWth 1DGF FWHM [log] SUM   #6
        flx_hdr     = 'S2G_FLS'#'S2G_FLS'
        lum_hdr     = 'S2G_LMS'#'S2G_LMF'
        flx_hdr_e   = 'S2G_FSE'#'S2G_FSE'
        lum_hdr_e1  = 'S2L_S'+str(error)+'L'# 'S2L_F'+str(error)+'L'
        lum_hdr_e2  = 'S2L_S'+str(error)+'H'# 'S2L_F'+str(error)+'H'

        lum_hdr_e11 = 'S2L_S1L'#'S2L_F1L'
        lum_hdr_e21 = 'S2L_S1H'#'S2L_F1H'
        lum_hdr_e12 = 'S2L_S2L'#'S2L_F2L'
        lum_hdr_e22 = 'S2L_S2H'#'S2L_F2H'
        lum_hdr_e13 = 'S2L_S3L'#'S2L_F3L'
        lum_hdr_e23 = 'S2L_S3H'#'S2L_F3H'

    elif method == 5 and log_lm == 'error':
        flx_hdr     = 'S2G_FLS'#'S2G_FLS'
        lum_hdr     = 'S2G_LMS'#'S2G_LMF'
        flx_hdr_e   = 'S2G_FSE'#'S2G_FSE'
        lum_hdr_e1  = 'S2M_S'+str(error)+'L'#'S2M_F'+str(error)+'L'
        lum_hdr_e2  = 'S2M_S'+str(error)+'L'#'S2M_F'+str(error)+'H'

        lum_hdr_e11 = 'S2M_S1L'#'S2M_F1L'
        lum_hdr_e21 = 'S2M_S1H'#'S2M_F1H'
        lum_hdr_e12 = 'S2M_S2L'#'S2M_F2L'
        lum_hdr_e22 = 'S2M_S2H'#'S2M_F2H'
        lum_hdr_e13 = 'S2M_S3L'#'S2M_F3L'
        lum_hdr_e23 = 'S2M_S3H'#'S2M_F3H'

    elif method == 5 and log_lm == 'value':
        flx_hdr     = 'S2G_FLS'#'S2G_FLS'
        lum_hdr     = 'S2G_LLS'#'S2G_LLF'
        flx_hdr_e   = 'S2G_FSE'#'S2G_FSE'
        lum_hdr_e1  = 'S2L_S' + str(error) + 'L'# 'S2L_F'+str(error)+'L'
        lum_hdr_e2  = 'S2L_S' + str(error) + 'H'# 'S2L_F'+str(error)+'H'

        lum_hdr_e11 = 'S2L_S1L'#'S2M_F1L'
        lum_hdr_e21 = 'S2L_S1H'#'S2M_F1H'
        lum_hdr_e12 = 'S2L_S2L'#'S2M_F2L'
        lum_hdr_e22 = 'S2L_S2H'#'S2M_F2H'
        lum_hdr_e13 = 'S2L_S3L'#'S2M_F3L'
        lum_hdr_e23 = 'S2L_S3H'#'S2M_F3H'

    elif method == 5 and log_lm == 'both':
        flx_hdr     = 'S2G_FLS'#'S2G_FLS'
        lum_hdr     = 'S2G_LLS'#'S2G_LLF'
        flx_hdr_e   = 'S2G_FSE'#'S2G_FSE'
        lum_hdr_e1  = 'S2M_S' + str(error) + 'L'#'S2M_F'+str(error)+'L'
        lum_hdr_e2  = 'S2M_S' + str(error) + 'H'#'S2M_F'+str(error)+'H'

        lum_hdr_e11 = 'S2M_S1L'#'S2M_F1L'
        lum_hdr_e21 = 'S2M_S1H'#'S2M_F1H'
        lum_hdr_e12 = 'S2M_S2L'#'S2M_F2L'
        lum_hdr_e22 = 'S2M_S2H'#'S2M_F2H'
        lum_hdr_e13 = 'S2M_S3L'#'S2M_F3L'
        lum_hdr_e23 = 'S2M_S3H'#'S2M_F3H'
    #DEFINES HEADERS TO READ FLX, LUM AND ERRS#

    ##########################TABLES########################
    ###########################STAT#########################
    SMPL = []
    N1_TOT         ,  N2_TOT        = [], []
    Z1_MED         ,  Z2_MED        = [], []
    V1_MED         ,  V2_MED        = [], []
    L1_AVG         ,  L2_AVG        = [], []
    L1_MED         ,  L2_MED        = [], []
    S1_AVG         ,  S2_AVG        = [], []
    S1_MED         ,  S2_MED        = [], []

    Z1_MED_E1      ,  Z2_MED_E1     = [], []
    Z1_MED_E2      ,  Z2_MED_E2     = [], []
    V1_MED_E1      ,  V2_MED_E1     = [], []
    V1_MED_E2      ,  V2_MED_E2     = [], []

    S1_AVG_E       ,  S2_AVG_E      = [], []
    S1_MED_E       ,  S2_MED_E      = [], []

    L1_AVG_E1      ,  L2_AVG_E1     = [], []
    L1_AVG_E2      ,  L2_AVG_E2     = [], []
    L1_MED_E1      ,  L2_MED_E1     = [], []
    L1_MED_E2      ,  L2_MED_E2     = [], []

    Z1_MED_E11     , Z1_MED_E21      = [], []
    Z1_MED_E12     , Z1_MED_E22      = [], []
    Z1_MED_E13     , Z1_MED_E23      = [], []
    V1_MED_E11     , V1_MED_E21      = [], []
    V1_MED_E12     , V1_MED_E22      = [], []
    V1_MED_E13     , V1_MED_E23      = [], []
    L1_AVG_E11     , L1_AVG_E21      = [], []
    L1_AVG_E12     , L1_AVG_E22      = [], []
    L1_AVG_E13     , L1_AVG_E23      = [], []
    L1_MED_E11     , L1_MED_E21      = [], []
    L1_MED_E12     , L1_MED_E22      = [], []
    L1_MED_E13     , L1_MED_E23      = [], []
    Z2_MED_E11     , Z2_MED_E21      = [], []
    Z2_MED_E12     , Z2_MED_E22      = [], []
    Z2_MED_E13     , Z2_MED_E23      = [], []
    V2_MED_E11     , V2_MED_E21      = [], []
    V2_MED_E12     , V2_MED_E22      = [], []
    V2_MED_E13     , V2_MED_E23      = [], []
    L2_AVG_E11     , L2_AVG_E21      = [], []
    L2_AVG_E12     , L2_AVG_E22      = [], []
    L2_AVG_E13     , L2_AVG_E23      = [], []
    L2_MED_E11     , L2_MED_E21      = [], []
    L2_MED_E12     , L2_MED_E22      = [], []
    L2_MED_E13     , L2_MED_E23      = [], []
    ###########################STAT#########################
    #################ARRAYS FOR TABLES FOR MCMC#############
    FA1_MSR_MCMC_OPT             = []
    FA1_MSR_MC_1_SG_L_MCMC_OPT,   FA1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    FA1_MSR_MC_2_SG_L_MCMC_OPT,   FA1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    FA1_MSR_MC_3_SG_L_MCMC_OPT,   FA1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    FA2_MSR_MCMC_OPT             = []
    FA2_MSR_MC_1_SG_L_MCMC_OPT,   FA2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    FA2_MSR_MC_2_SG_L_MCMC_OPT,   FA2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    FA2_MSR_MC_3_SG_L_MCMC_OPT,   FA2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    FA12_MSR_MCMC_OPT            = []
    FA12_MSR_MC_1_SG_L_MCMC_OPT,  FA12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    FA12_MSR_MC_2_SG_L_MCMC_OPT,  FA12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    FA12_MSR_MC_3_SG_L_MCMC_OPT,  FA12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    FM1_MSR_MCMC_OPT             = []
    FM1_MSR_MC_1_SG_L_MCMC_OPT,   FM1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    FM1_MSR_MC_2_SG_L_MCMC_OPT,   FM1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    FM1_MSR_MC_3_SG_L_MCMC_OPT,   FM1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    FM2_MSR_MCMC_OPT             = []
    FM2_MSR_MC_1_SG_L_MCMC_OPT,   FM2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    FM2_MSR_MC_2_SG_L_MCMC_OPT,   FM2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    FM2_MSR_MC_3_SG_L_MCMC_OPT,   FM2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    FM12_MSR_MCMC_OPT            = []
    FM12_MSR_MC_1_SG_L_MCMC_OPT,  FM12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    FM12_MSR_MC_2_SG_L_MCMC_OPT,  FM12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    FM12_MSR_MC_3_SG_L_MCMC_OPT,  FM12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LA1_MSR_MCMC_OPT             = []
    LA1_MSR_MC_1_SG_L_MCMC_OPT,   LA1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    LA1_MSR_MC_2_SG_L_MCMC_OPT,   LA1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    LA1_MSR_MC_3_SG_L_MCMC_OPT,   LA1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    LA2_MSR_MCMC_OPT             = []
    LA2_MSR_MC_1_SG_L_MCMC_OPT,   LA2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    LA2_MSR_MC_2_SG_L_MCMC_OPT,   LA2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    LA2_MSR_MC_3_SG_L_MCMC_OPT,   LA2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    LA12_MSR_MCMC_OPT            = []
    LA12_MSR_MC_1_SG_L_MCMC_OPT,  LA12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LA12_MSR_MC_2_SG_L_MCMC_OPT,  LA12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LA12_MSR_MC_3_SG_L_MCMC_OPT,  LA12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LM1_MSR_MCMC_OPT             = []
    LM1_MSR_MC_1_SG_L_MCMC_OPT,   LM1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    LM1_MSR_MC_2_SG_L_MCMC_OPT,   LM1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    LM1_MSR_MC_3_SG_L_MCMC_OPT,   LM1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    LM2_MSR_MCMC_OPT             = []
    LM2_MSR_MC_1_SG_L_MCMC_OPT,   LM2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    LM2_MSR_MC_2_SG_L_MCMC_OPT,   LM2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    LM2_MSR_MC_3_SG_L_MCMC_OPT,   LM2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    LM12_MSR_MCMC_OPT            = []
    LM12_MSR_MC_1_SG_L_MCMC_OPT,  LM12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LM12_MSR_MC_2_SG_L_MCMC_OPT,  LM12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LM12_MSR_MC_3_SG_L_MCMC_OPT,  LM12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LLA1_MSR_MCMC_OPT            = []
    LLA1_MSR_MC_1_SG_L_MCMC_OPT,  LLA1_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LLA1_MSR_MC_2_SG_L_MCMC_OPT,  LLA1_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LLA1_MSR_MC_3_SG_L_MCMC_OPT,  LLA1_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LLA2_MSR_MCMC_OPT            = []
    LLA2_MSR_MC_1_SG_L_MCMC_OPT,  LLA2_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LLA2_MSR_MC_2_SG_L_MCMC_OPT,  LLA2_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LLA2_MSR_MC_3_SG_L_MCMC_OPT,  LLA2_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LLA12_MSR_MCMC_OPT           = []
    LLA12_MSR_MC_1_SG_L_MCMC_OPT, LLA12_MSR_MC_1_SG_H_MCMC_OPT = [], []
    LLA12_MSR_MC_2_SG_L_MCMC_OPT, LLA12_MSR_MC_2_SG_H_MCMC_OPT = [], []
    LLA12_MSR_MC_3_SG_L_MCMC_OPT, LLA12_MSR_MC_3_SG_H_MCMC_OPT = [], []
    LLM1_MSR_MCMC_OPT            = []
    LLM1_MSR_MC_1_SG_L_MCMC_OPT,  LLM1_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LLM1_MSR_MC_2_SG_L_MCMC_OPT,  LLM1_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LLM1_MSR_MC_3_SG_L_MCMC_OPT,  LLM1_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LLM2_MSR_MCMC_OPT            = []
    LLM2_MSR_MC_1_SG_L_MCMC_OPT,  LLM2_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    LLM2_MSR_MC_2_SG_L_MCMC_OPT,  LLM2_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    LLM2_MSR_MC_3_SG_L_MCMC_OPT,  LLM2_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    LLM12_MSR_MCMC_OPT           = []
    LLM12_MSR_MC_1_SG_L_MCMC_OPT, LLM12_MSR_MC_1_SG_H_MCMC_OPT = [], []
    LLM12_MSR_MC_2_SG_L_MCMC_OPT, LLM12_MSR_MC_2_SG_H_MCMC_OPT = [], []
    LLM12_MSR_MC_3_SG_L_MCMC_OPT, LLM12_MSR_MC_3_SG_H_MCMC_OPT = [], []
    VA1_MSR_MCMC_OPT             = []
    VA1_MSR_MC_1_SG_L_MCMC_OPT,   VA1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    VA1_MSR_MC_2_SG_L_MCMC_OPT,   VA1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    VA1_MSR_MC_3_SG_L_MCMC_OPT,   VA1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    VA2_MSR_MCMC_OPT             = []
    VA2_MSR_MC_1_SG_L_MCMC_OPT,   VA2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    VA2_MSR_MC_2_SG_L_MCMC_OPT,   VA2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    VA2_MSR_MC_3_SG_L_MCMC_OPT,   VA2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    VA12_MSR_MCMC_OPT            = []
    VA12_MSR_MC_1_SG_L_MCMC_OPT,  VA12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    VA12_MSR_MC_2_SG_L_MCMC_OPT,  VA12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    VA12_MSR_MC_3_SG_L_MCMC_OPT,  VA12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    VM1_MSR_MCMC_OPT             = []
    VM1_MSR_MC_1_SG_L_MCMC_OPT,   VM1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    VM1_MSR_MC_2_SG_L_MCMC_OPT,   VM1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    VM1_MSR_MC_3_SG_L_MCMC_OPT,   VM1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    VM2_MSR_MCMC_OPT             = []
    VM2_MSR_MC_1_SG_L_MCMC_OPT,   VM2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    VM2_MSR_MC_2_SG_L_MCMC_OPT,   VM2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    VM2_MSR_MC_3_SG_L_MCMC_OPT,   VM2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    VM12_MSR_MCMC_OPT            = []
    VM12_MSR_MC_1_SG_L_MCMC_OPT,  VM12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    VM12_MSR_MC_2_SG_L_MCMC_OPT,  VM12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    VM12_MSR_MC_3_SG_L_MCMC_OPT,  VM12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    ZA1_MSR_MCMC_OPT             = []
    ZA1_MSR_MC_1_SG_L_MCMC_OPT,   ZA1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    ZA1_MSR_MC_2_SG_L_MCMC_OPT,   ZA1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    ZA1_MSR_MC_3_SG_L_MCMC_OPT,   ZA1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    ZA2_MSR_MCMC_OPT             = []
    ZA2_MSR_MC_1_SG_L_MCMC_OPT,   ZA2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    ZA2_MSR_MC_2_SG_L_MCMC_OPT,   ZA2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    ZA2_MSR_MC_3_SG_L_MCMC_OPT,   ZA2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    ZA12_MSR_MCMC_OPT            = []
    ZA12_MSR_MC_1_SG_L_MCMC_OPT,  ZA12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    ZA12_MSR_MC_2_SG_L_MCMC_OPT,  ZA12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    ZA12_MSR_MC_3_SG_L_MCMC_OPT,  ZA12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    ZM1_MSR_MCMC_OPT             = []
    ZM1_MSR_MC_1_SG_L_MCMC_OPT,   ZM1_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    ZM1_MSR_MC_2_SG_L_MCMC_OPT,   ZM1_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    ZM1_MSR_MC_3_SG_L_MCMC_OPT,   ZM1_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    ZM2_MSR_MCMC_OPT             = []
    ZM2_MSR_MC_1_SG_L_MCMC_OPT,   ZM2_MSR_MC_1_SG_H_MCMC_OPT   = [], []
    ZM2_MSR_MC_2_SG_L_MCMC_OPT,   ZM2_MSR_MC_2_SG_H_MCMC_OPT   = [], []
    ZM2_MSR_MC_3_SG_L_MCMC_OPT,   ZM2_MSR_MC_3_SG_H_MCMC_OPT   = [], []
    ZM12_MSR_MCMC_OPT            = []
    ZM12_MSR_MC_1_SG_L_MCMC_OPT,  ZM12_MSR_MC_1_SG_H_MCMC_OPT  = [], []
    ZM12_MSR_MC_2_SG_L_MCMC_OPT,  ZM12_MSR_MC_2_SG_H_MCMC_OPT  = [], []
    ZM12_MSR_MC_3_SG_L_MCMC_OPT,  ZM12_MSR_MC_3_SG_H_MCMC_OPT  = [], []
    #################ARRAYS FOR TABLES FOR MCMC#############

    ##############ARRAYS FOR TABLES FOR MCMC PLOT###########
    Z1_AVG_E1_MC_OPT    , Z1_AVG_E2_MC_OPT    = [], []
    Z2_AVG_E1_MC_OPT    , Z2_AVG_E2_MC_OPT    = [], []
    Z1_MED_E1_MC_OPT    , Z1_MED_E2_MC_OPT    = [], []
    Z2_MED_E1_MC_OPT    , Z2_MED_E2_MC_OPT    = [], []
    V1_AVG_E1_MC_OPT    , V1_AVG_E2_MC_OPT    = [], []
    V2_AVG_E1_MC_OPT    , V2_AVG_E2_MC_OPT    = [], []
    V1_MED_E1_MC_OPT    , V1_MED_E2_MC_OPT    = [], []
    V2_MED_E1_MC_OPT    , V2_MED_E2_MC_OPT    = [], []
    S1_AVG_E1_MC_OPT    , S1_AVG_E2_MC_OPT    = [], []
    S2_AVG_E1_MC_OPT    , S2_AVG_E2_MC_OPT    = [], []
    S1_MED_E1_MC_OPT    , S1_MED_E2_MC_OPT    = [], []
    S2_MED_E1_MC_OPT    , S2_MED_E2_MC_OPT    = [], []
    L1_AVG_E1_MC_1_OPT  , L1_AVG_E2_MC_1_OPT  = [], []
    L1_MED_E1_MC_1_OPT  , L1_MED_E2_MC_1_OPT  = [], []
    L1_AVG_E1_MC_2_OPT  , L1_AVG_E2_MC_2_OPT  = [], []
    L1_MED_E1_MC_2_OPT  , L1_MED_E2_MC_2_OPT  = [], []
    L2_AVG_E1_MC_1_OPT  , L2_AVG_E2_MC_1_OPT  = [], []
    L2_MED_E1_MC_1_OPT  , L2_MED_E2_MC_1_OPT  = [], []
    L2_AVG_E1_MC_2_OPT  , L2_AVG_E2_MC_2_OPT  = [], []
    L2_MED_E1_MC_2_OPT  , L2_MED_E2_MC_2_OPT  = [], []
    S12_AVG_E1_MC_OPT   , S12_AVG_E2_MC_OPT   = [], []
    S12_MED_E1_MC_OPT   , S12_MED_E2_MC_OPT   = [], []
    L12_AVG_E1_MC_1_OPT , L12_AVG_E2_MC_1_OPT = [], []
    L12_MED_E1_MC_1_OPT , L12_MED_E2_MC_1_OPT = [], []
    L12_AVG_E1_MC_2_OPT , L12_AVG_E2_MC_2_OPT = [], []
    L12_MED_E1_MC_2_OPT , L12_MED_E2_MC_2_OPT = [], []
    ##############ARRAYS FOR TABLES FOR MCMC PLOT###########

    ##########################TABLES########################
    for element in sbsmn_mcmc:

        prefix_line1    = line1 + '-'
        prefix_line2    = line2 + '-'
        stk_hme_dir_1   = home + 'Stack_Results-'+ line1 +'-3D/'
        img_dir_res_1   = stk_hme_dir_1 + 'IMAGES/'  
        stp_dir_res_1   = stk_hme_dir_1 + 'STAMPS/'  + str(spc_wdt_dir) +'/' 
        tbl_dir_res_1   = stk_hme_dir_1 + 'TABLES/'  + str(spc_wdt_dir) +'/'
        mcm_dir_tbl_1   = tbl_dir_res_1 + 'MCMC/'
        plt_dir_tbl_1   = tbl_dir_res_1 + 'PLOTS/'  
        plt_dir_res_1   = stk_hme_dir_1 + 'PLOTS/'   + str(spc_wdt_dir) +'/'
        mcm_dir_plt_1   = plt_dir_res_1 + 'MCMC/'    
        res_dir_plt_1   = plt_dir_res_1 + 'RESULTS/'    
        stk_dir_res_1   = stk_hme_dir_1 + 'STACKS/'  + str(spc_wdt_dir) +'/'

        if line1 == '13CO':
            restframe_frequency_1      =   110.20137E9           
        elif line1 == '12CO':
            restframe_frequency_1      =   115.271208E9
        elif line1 == '18CO':
            restframe_frequency_1      =   109.78217340E9


        stk_hme_dir_2   = home + 'Stack_Results-'+ line2 +'-3D/'
        img_dir_res_2   = stk_hme_dir_2 + 'IMAGES/' + str(spc_wdt_dir) +'/' 
        stp_dir_res_2   = stk_hme_dir_2 + 'STAMPS/' + str(spc_wdt_dir) +'/' 
        tbl_dir_res_2   = stk_hme_dir_2 + 'TABLES/' + str(spc_wdt_dir) +'/'
        mcm_dir_tbl_2   = tbl_dir_res_2 + 'MCMC/'
        plt_dir_tbl_2   = tbl_dir_res_2 + 'PLOTS/'  
        plt_dir_res_2   = stk_hme_dir_2 + 'PLOTS/'  + str(spc_wdt_dir) +'/'
        mcm_dir_plt_2   = plt_dir_res_2 + 'MCMC/'    
        res_dir_plt_2   = plt_dir_res_2 + 'RESULTS/'
        stk_dir_res_2   = stk_hme_dir_2 + 'STACKS/' + str(spc_wdt_dir) +'/'
        if line2 == '13CO':
            restframe_frequency_2      =   110.20137E9           
        elif line2 == '12CO':
            restframe_frequency_2      =   115.271208E9
        elif line2 == '18CO':
            restframe_frequency_2      =   109.78217340E9

        DIR_RES_1     = [
                            stk_hme_dir_1,img_dir_res_1,stp_dir_res_1,
                            tbl_dir_res_1,mcm_dir_tbl_1,plt_dir_tbl_1,
                            plt_dir_res_1,mcm_dir_plt_1,res_dir_plt_1,
                            stk_dir_res_1
                        ]
        DIR_RES_2     = [
                            stk_hme_dir_2,img_dir_res_2,stp_dir_res_2,
                            tbl_dir_res_2,mcm_dir_tbl_2,plt_dir_tbl_2,
                            plt_dir_res_2,mcm_dir_plt_2,res_dir_plt_2,
                            stk_dir_res_2
                        ]               
        cat_ipt_tbl   =   cat_dir + 'CII_Sources_HATLAS-' + line + '-' + str(sbsms_mcmc) + '-' +str(element) 
        print 
        print cat_ipt_tbl
        print
        Check_directories(cat_ipt_tbl,cat_parent,DIR_RES=DIR_RES_1)
        Check_directories(cat_ipt_tbl,cat_parent,DIR_RES=DIR_RES_2)        
        spec_file_plt_red_1_avg = stp_dir_res_1 + prefix_line1 + 'CII_HATLAS-' + str(sbsms_mcmc) + '-' +str(element) +'-stk-'+func1 + '-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_dta_ms-2DC-sum.fits'
        spec_file_plt_red_1_med = stp_dir_res_1 + prefix_line1 + 'CII_HATLAS-' + str(sbsms_mcmc) + '-' +str(element) +'-stk-'+func2 + '-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_dta_ms-2DC-sum.fits'
        
        spec_file_plt_red_2_avg = stp_dir_res_2 + prefix_line2 + 'CII_HATLAS-' + str(sbsms_mcmc) + '-' +str(element) +'-stk-'+func1 + '-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_dta_ms-2DC-sum.fits'
        spec_file_plt_red_2_med = stp_dir_res_2 + prefix_line2 + 'CII_HATLAS-' + str(sbsms_mcmc) + '-' +str(element) +'-stk-'+func2 + '-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_dta_ms-2DC-sum.fits'

        cat_ipt_tbl_1  = cat_dir + 'CII_Sources_HATLAS-' + line1 + '-' + str(sbsms_mcmc) + '-' + str(element)
        cat_tbl_1      = cat_ipt_tbl_1 + tbl_ext_ipt
        Cat_Ipt_Tbl_1  = Table_Read(cat_tbl_1,tbl_format_ipt)
        z_1            = Cat_Ipt_Tbl_1[8]
        v_1            = Cat_Ipt_Tbl_1[Splt_CNm]

        cat_ipt_tbl_2  = cat_dir + 'CII_Sources_HATLAS-' + line2 + '-' + str(sbsms_mcmc) + '-' + str(element)
        cat_tbl_2      = cat_ipt_tbl_2 + tbl_ext_ipt
        Cat_Ipt_Tbl_2  = Table_Read(cat_tbl_2,tbl_format_ipt)
        z_2            = Cat_Ipt_Tbl_2[8]
        v_2            = Cat_Ipt_Tbl_2[Splt_CNm]

        print
        print ('1: ',spec_file_plt_red_1_avg)
        print ('1: ',spec_file_plt_red_1_med)
        print ('2: ',spec_file_plt_red_2_avg)
        print ('2: ',spec_file_plt_red_2_med)
        print
        print ('Catalogue : ',cat_tbl_1)
        print ('Catalogue : ',cat_tbl_2)
        print

        ntt_hdr     = 'STK_NUM'
        rds_hdr     = 'STZ_MED'
        rds_hdr_e1  = 'STZ_'+str(error)+'SL'
        rds_hdr_e2  = 'STZ_'+str(error)+'SH'

        rds_hdr_e11 = 'STZ_1SL'
        rds_hdr_e21 = 'STZ_1SH'
        rds_hdr_e12 = 'STZ_2SL'
        rds_hdr_e22 = 'STZ_2SH'
        rds_hdr_e13 = 'STZ_3SL'
        rds_hdr_e23 = 'STZ_3SH'


        ntt_hdr     = 'STK_NUM'
        vrx_hdr     = 'STS_MED'
        vrx_hdr_e1  = 'STS_'+str(error)+'SL'
        vrx_hdr_e2  = 'STS_'+str(error)+'SH'

        vrx_hdr_e11 = 'STS_1SL'
        vrx_hdr_e21 = 'STS_1SH'
        vrx_hdr_e12 = 'STS_2SL'
        vrx_hdr_e22 = 'STS_2SH'
        vrx_hdr_e13 = 'STS_3SL'
        vrx_hdr_e23 = 'STS_3SH'

        SMPL.append(str(element))

        N1_TOT.append(Header_Get(spec_file_plt_red_1_avg,ntt_hdr))
        Z1_MED.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr))
        V1_MED.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr))

        S1_AVG.append(Header_Get(spec_file_plt_red_1_avg,flx_hdr))
        S1_MED.append(Header_Get(spec_file_plt_red_1_med,flx_hdr))

        L1_AVG.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr))
        L1_MED.append(Header_Get(spec_file_plt_red_1_med,lum_hdr))

        Z1_MED_E1.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e1))
        Z1_MED_E2.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e2))
        V1_MED_E1.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e1))
        V1_MED_E2.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e2))

        S1_AVG_E.append(Header_Get(spec_file_plt_red_1_avg,flx_hdr_e))
        S1_MED_E.append(Header_Get(spec_file_plt_red_1_med,flx_hdr_e))

        L1_AVG_E1.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e1))
        L1_AVG_E2.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e2))
        L1_MED_E1.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e1))
        L1_MED_E2.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e2))

        N2_TOT.append(Header_Get(spec_file_plt_red_2_avg,ntt_hdr))
        Z2_MED.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr))
        V2_MED.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr))

        S2_AVG.append(Header_Get(spec_file_plt_red_2_avg,flx_hdr))
        S2_MED.append(Header_Get(spec_file_plt_red_2_med,flx_hdr))

        L2_AVG.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr))
        L2_MED.append(Header_Get(spec_file_plt_red_2_med,lum_hdr))

        Z2_MED_E1.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e1))
        Z2_MED_E2.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e2))
        V2_MED_E1.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e1))
        V2_MED_E2.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e2))

        S2_AVG_E.append(Header_Get(spec_file_plt_red_2_avg,flx_hdr_e))
        S2_MED_E.append(Header_Get(spec_file_plt_red_2_med,flx_hdr_e))

        L2_AVG_E1.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e1))
        L2_AVG_E2.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e2))
        L2_MED_E1.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e1))
        L2_MED_E2.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e2))

        Z1_MED_E11.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e11))
        Z1_MED_E21.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e21))
        Z1_MED_E12.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e12))
        Z1_MED_E22.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e22))
        Z1_MED_E13.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e13))
        Z1_MED_E23.append(Header_Get(spec_file_plt_red_1_avg,rds_hdr_e23))

        V1_MED_E11.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e11))
        V1_MED_E21.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e21))
        V1_MED_E12.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e12))
        V1_MED_E22.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e22))
        V1_MED_E13.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e13))
        V1_MED_E23.append(Header_Get(spec_file_plt_red_1_avg,vrx_hdr_e23))

        L1_AVG_E11.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e11))
        L1_AVG_E21.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e21))
        L1_AVG_E12.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e12))
        L1_AVG_E22.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e22))
        L1_AVG_E13.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e13))
        L1_AVG_E23.append(Header_Get(spec_file_plt_red_1_avg,lum_hdr_e23))

        L1_MED_E11.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e11))
        L1_MED_E21.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e21))
        L1_MED_E12.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e12))
        L1_MED_E22.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e22))
        L1_MED_E13.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e13))
        L1_MED_E23.append(Header_Get(spec_file_plt_red_1_med,lum_hdr_e23))

        Z2_MED_E11.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e11))
        Z2_MED_E21.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e21))
        Z2_MED_E12.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e12))
        Z2_MED_E22.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e22))
        Z2_MED_E13.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e13))
        Z2_MED_E23.append(Header_Get(spec_file_plt_red_2_avg,rds_hdr_e23))

        V2_MED_E11.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e11))
        V2_MED_E21.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e21))
        V2_MED_E12.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e12))
        V2_MED_E22.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e22))
        V2_MED_E13.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e13))
        V2_MED_E23.append(Header_Get(spec_file_plt_red_2_avg,vrx_hdr_e23))

        L2_AVG_E11.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e11))
        L2_AVG_E21.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e21))
        L2_AVG_E12.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e12))
        L2_AVG_E22.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e22))
        L2_AVG_E13.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e13))
        L2_AVG_E23.append(Header_Get(spec_file_plt_red_2_avg,lum_hdr_e23))

        L2_MED_E11.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e11))
        L2_MED_E21.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e21))
        L2_MED_E12.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e12))
        L2_MED_E22.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e22))
        L2_MED_E13.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e13))
        L2_MED_E23.append(Header_Get(spec_file_plt_red_2_med,lum_hdr_e23))

        #####################TABLE-01#####################
        mct01 = aptbl.Table()
        mct01['SMPL'] = SMPL

        mct01['N1_TOT']            = N1_TOT
        mct01['Z1_MED']            = Z1_MED
        mct01['Z1_MED_E1SGL']      = Z1_MED_E11
        mct01['Z1_MED_E1SGH']      = Z1_MED_E21
        mct01['Z1_MED_E2SGL']      = Z1_MED_E12
        mct01['Z1_MED_E2SGH']      = Z1_MED_E22
        mct01['Z1_MED_E3SGL']      = Z1_MED_E13
        mct01['Z1_MED_E3SGH']      = Z1_MED_E23

        mct01['N2_TOT']            = N2_TOT
        mct01['Z2_MED']            = Z2_MED
        mct01['Z2_MED_E1SGL']      = Z2_MED_E11
        mct01['Z2_MED_E1SGH']      = Z2_MED_E21
        mct01['Z2_MED_E2SGL']      = Z2_MED_E12
        mct01['Z2_MED_E2SGH']      = Z2_MED_E22
        mct01['Z2_MED_E3SGL']      = Z2_MED_E13
        mct01['Z2_MED_E3SGH']      = Z2_MED_E23

        TABLESTATNAME_1_1_1  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_1_1_2  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        TABLESTATNAME_1_2_1  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_1_2_2  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'      

        mct01.write(TABLESTATNAME_1_1_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct01.write(TABLESTATNAME_1_1_2, format=tbl_format_opt, overwrite = True)

        mct01.write(TABLESTATNAME_1_2_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct01.write(TABLESTATNAME_1_2_2, format=tbl_format_opt, overwrite = True)
        #####################TABLE-01#####################

        #####################TABLE-02#####################
        mct02 = aptbl.Table()
        mct02['SMPL'] = SMPL

        mct02['N1_TOT']            = N1_TOT
        mct02['S1_AVG']            = S1_AVG
        mct02['S1_AVG_E']          = S1_AVG_E
        mct02['S1_MED']            = S1_MED
        mct02['S1_MED_E']          = S1_MED_E

        mct02['N2_TOT']            = N2_TOT
        mct02['S2_AVG']            = S2_AVG
        mct02['S2_AVG_E']          = S2_AVG_E
        mct02['S2_MED']            = S2_MED
        mct02['S2_MED_E']          = S2_MED_E
        TABLESTATNAME_2_1_1  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_2_1_2  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        TABLESTATNAME_2_2_1  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_2_2_2  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct02.write(TABLESTATNAME_2_1_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct02.write(TABLESTATNAME_2_1_2, format=tbl_format_opt, overwrite = True)

        mct02.write(TABLESTATNAME_2_2_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct02.write(TABLESTATNAME_2_2_2, format=tbl_format_opt, overwrite = True)
        #####################TABLE-02#####################

        #####################TABLE-03#####################
        mct03 = aptbl.Table()
        mct03['SMPL'] = SMPL

        mct03['N1_TOT']            = N1_TOT
        mct03['L1_AVG']            = L1_AVG
        mct03['L1_AVG_E1SGL']      = L1_AVG_E11
        mct03['L1_AVG_E1SGH']      = L1_AVG_E21
        mct03['L1_AVG_E2SGL']      = L1_AVG_E12
        mct03['L1_AVG_E2SGH']      = L1_AVG_E22
        mct03['L1_AVG_E3SGL']      = L1_AVG_E13
        mct03['L1_AVG_E3SGH']      = L1_AVG_E23

        mct03['L1_MED']            = L1_MED
        mct03['L1_MED_E1SGL']      = L1_MED_E11
        mct03['L1_MED_E1SGH']      = L1_MED_E21
        mct03['L1_MED_E2SGL']      = L1_MED_E12
        mct03['L1_MED_E2SGH']      = L1_MED_E22
        mct03['L1_MED_E3SGL']      = L1_MED_E13
        mct03['L1_MED_E3SGH']      = L1_MED_E23

        mct03['N2_TOT']            = N2_TOT
        mct03['L2_AVG']            = L2_AVG
        mct03['L2_AVG_E1SGL']      = L2_AVG_E11
        mct03['L2_AVG_E1SGH']      = L2_AVG_E21
        mct03['L2_AVG_E2SGL']      = L2_AVG_E12
        mct03['L2_AVG_E2SGH']      = L2_AVG_E22
        mct03['L2_AVG_E3SGL']      = L2_AVG_E13
        mct03['L2_AVG_E3SGH']      = L2_AVG_E23

        mct03['L2_MED']            = L2_MED
        mct03['L2_MED_E1SGL']      = L2_MED_E11
        mct03['L2_MED_E1SGH']      = L2_MED_E21
        mct03['L2_MED_E2SGL']      = L2_MED_E12
        mct03['L2_MED_E2SGH']      = L2_MED_E22
        mct03['L2_MED_E3SGL']      = L2_MED_E13
        mct03['L2_MED_E3SGH']      = L2_MED_E23

        TABLESTATNAME_3_1_1  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_3_1_2  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        TABLESTATNAME_3_2_1  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_3_2_2  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct03.write(TABLESTATNAME_3_1_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct03.write(TABLESTATNAME_3_1_2, format=tbl_format_opt, overwrite = True)

        mct03.write(TABLESTATNAME_3_2_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct03.write(TABLESTATNAME_3_2_2, format=tbl_format_opt, overwrite = True)
        #####################TABLE-03#####################

        #####################TABLE-04#####################
        mct04 = aptbl.Table()
        mct04['SMPL'] = SMPL

        mct04['N1_TOT']            = N1_TOT
        mct04['V1_MED']            = V1_MED
        mct04['V1_MED_E1SGL']      = V1_MED_E11
        mct04['V1_MED_E1SGH']      = V1_MED_E21
        mct04['V1_MED_E2SGL']      = V1_MED_E12
        mct04['V1_MED_E2SGH']      = V1_MED_E22
        mct04['V1_MED_E3SGL']      = V1_MED_E13
        mct04['V1_MED_E3SGH']      = V1_MED_E23

        mct04['N2_TOT']            = N2_TOT
        mct04['V2_MED']            = V2_MED
        mct04['V2_MED_E1SGL']      = V2_MED_E11
        mct04['V2_MED_E1SGH']      = V2_MED_E21
        mct04['V2_MED_E2SGL']      = V2_MED_E12
        mct04['V2_MED_E2SGH']      = V2_MED_E22
        mct04['V2_MED_E3SGL']      = V2_MED_E13
        mct04['V2_MED_E3SGH']      = V2_MED_E23
        TABLESTATNAME_4_1_1  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_4_1_2  = mcm_dir_tbl_1 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        TABLESTATNAME_4_2_1  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_4_2_2  = mcm_dir_tbl_2 + '/CII_HATLAS-' + str(sbsms_mcmc) + '-' + 'MS-'+str(method)+'-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct04.write(TABLESTATNAME_4_1_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct04.write(TABLESTATNAME_4_1_2, format=tbl_format_opt, overwrite = True)

        mct04.write(TABLESTATNAME_4_2_1, format='ascii.fixed_width_two_line', overwrite = True)
        mct04.write(TABLESTATNAME_4_2_2, format=tbl_format_opt, overwrite = True)
        ####################TABLE-04#####################

        ################################################MCMC-COMP#################################################
        MCMC_CDF_VLS_PLT_OPT = MCMC_Confidence_Interval(iterations_mc,line1,line2,method,error,
                                                    z_1,v_1,z_2,v_2,
                                                    var_smpls_mcmc = sbsms_mcmc,
                                                    nmb_smpls_mcmc = element,
                                                    z_lne1_mcmc = z_1,z_lne2_mcmc = z_2,
                                                    v_lne1_mcmc = v_1,v_lne2_mcmc = v_2,
                                                    spec_file_plt_lne1_avg_mcmc = spec_file_plt_red_1_avg,
                                                    spec_file_plt_lne1_med_mcmc = spec_file_plt_red_1_med,
                                                    spec_file_plt_lne2_avg_mcmc = spec_file_plt_red_2_avg,
                                                    spec_file_plt_lne2_med_mcmc = spec_file_plt_red_2_med,
                                                    flx_lne1_avg_hdr   = flx_hdr  ,flx_lne1_med_hdr   = flx_hdr,
                                                    flx_lne1_avg_hdr_e = flx_hdr_e,flx_lne1_med_hdr_e = flx_hdr_e,
                                                    flx_lne2_avg_hdr   = flx_hdr  ,flx_lne2_med_hdr   = flx_hdr,
                                                    flx_lne2_avg_hdr_e = flx_hdr_e,flx_lne2_med_hdr_e = flx_hdr_e,
                                                    flx_hdr     = flx_hdr,
                                                    lum_hdr     = lum_hdr,
                                                    flx_hdr_e   = flx_hdr_e,
                                                    lum_hdr_e1  = lum_hdr_e1,
                                                    lum_hdr_e2  = lum_hdr_e2,
                                                    lum_hdr_e11 = lum_hdr_e11,
                                                    lum_hdr_e21 = lum_hdr_e21,
                                                    lum_hdr_e12 = lum_hdr_e12,
                                                    lum_hdr_e22 = lum_hdr_e22,
                                                    lum_hdr_e13 = lum_hdr_e13,
                                                    lum_hdr_e23 = lum_hdr_e23,
                                                    plot_dist_hist = plot_dist_hist,
                                                    dest_dir_plt       = mcm_dir_plt_1,**kwargs)

        ################################################MCMC-COMP#################################################

        MCMC_CDF_VLS_OPT = MCMC_CDF_VLS_PLT_OPT[0]
        MCMC_CDF_PLT_OPT = MCMC_CDF_VLS_PLT_OPT[1]
        FA1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[0][0])
        FA1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[1][0])
        FA1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[2][0])
        FA1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[3][0])
        FA1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[4][0])
        FA1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[5][0])
        FA1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[6][0])
        FA2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[7][0])
        FA2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[8][0])
        FA2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[9][0])
        FA2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[10][0])
        FA2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[11][0])
        FA2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[12][0])
        FA2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[13][0])
        FA12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[14][0])
        FA12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[15][0])
        FA12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[16][0])
        FA12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[17][0])
        FA12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[18][0])
        FA12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[19][0])
        FA12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[20][0])
        FM1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[21][0])
        FM1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[22][0])
        FM1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[23][0])
        FM1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[24][0])
        FM1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[25][0])
        FM1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[26][0])
        FM1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[27][0])
        FM2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[28][0])
        FM2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[29][0])
        FM2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[30][0])
        FM2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[31][0])
        FM2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[32][0])
        FM2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[33][0])
        FM2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[34][0])
        FM12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[35][0])
        FM12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[36][0])
        FM12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[37][0])
        FM12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[38][0])
        FM12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[39][0])
        FM12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[40][0])
        FM12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[41][0])
        LA1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[42][0])
        LA1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[43][0])
        LA1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[44][0])
        LA1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[45][0])
        LA1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[46][0])
        LA1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[47][0])
        LA1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[48][0])
        LA2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[49][0])
        LA2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[50][0])
        LA2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[51][0])
        LA2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[52][0])
        LA2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[53][0])
        LA2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[54][0])
        LA2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[55][0])
        LA12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[56][0])
        LA12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[57][0])
        LA12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[58][0])
        LA12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[59][0])
        LA12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[60][0])
        LA12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[61][0])
        LA12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[62][0])
        LM1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[63][0])
        LM1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[64][0])
        LM1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[65][0])
        LM1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[66][0])
        LM1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[67][0])
        LM1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[68][0])
        LM1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[69][0])
        LM2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[70][0])
        LM2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[71][0])
        LM2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[72][0])
        LM2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[73][0])
        LM2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[74][0])
        LM2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[75][0])
        LM2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[76][0])
        LM12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[77][0])
        LM12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[78][0])
        LM12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[79][0])
        LM12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[80][0])
        LM12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[81][0])
        LM12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[82][0])
        LM12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[83][0])
        LLA1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[84][0])
        LLA1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[85][0])
        LLA1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[86][0])
        LLA1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[87][0])
        LLA1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[88][0])
        LLA1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[89][0])
        LLA1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[90][0])
        LLA2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[91][0])
        LLA2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[92][0])
        LLA2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[93][0])
        LLA2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[94][0])
        LLA2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[95][0])
        LLA2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[96][0])
        LLA2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[97][0])
        LLA12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[98][0])
        LLA12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[99][0])
        LLA12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[100][0])
        LLA12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[101][0])
        LLA12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[102][0])
        LLA12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[103][0])
        LLA12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[104][0])
        LLM1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[105][0])
        LLM1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[106][0])
        LLM1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[107][0])
        LLM1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[108][0])
        LLM1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[109][0])
        LLM1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[110][0])
        LLM1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[111][0])
        LLM2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[112][0])
        LLM2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[113][0])
        LLM2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[114][0])
        LLM2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[115][0])
        LLM2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[116][0])
        LLM2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[117][0])
        LLM2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[118][0])
        LLM12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[119][0])
        LLM12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[120][0])
        LLM12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[121][0])
        LLM12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[122][0])
        LLM12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[123][0])
        LLM12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[124][0])
        LLM12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[125][0])
        VA1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[126][0])
        VA1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[127][0])
        VA1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[128][0])
        VA1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[129][0])
        VA1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[130][0])
        VA1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[131][0])
        VA1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[132][0])
        VA2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[133][0])
        VA2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[134][0])
        VA2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[135][0])
        VA2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[136][0])
        VA2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[137][0])
        VA2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[138][0])
        VA2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[139][0])
        VA12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[140][0])
        VA12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[141][0])
        VA12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[142][0])
        VA12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[143][0])
        VA12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[144][0])
        VA12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[145][0])
        VA12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[146][0])
        VM1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[147][0])
        VM1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[148][0])
        VM1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[149][0])
        VM1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[150][0])
        VM1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[151][0])
        VM1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[152][0])
        VM1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[153][0])
        VM2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[154][0])
        VM2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[155][0])
        VM2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[156][0])
        VM2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[157][0])
        VM2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[158][0])
        VM2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[159][0])
        VM2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[160][0])
        VM12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[161][0])
        VM12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[162][0])
        VM12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[163][0])
        VM12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[164][0])
        VM12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[165][0])
        VM12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[166][0])
        VM12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[167][0])
        ZA1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[168][0])
        ZA1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[169][0])
        ZA1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[170][0])
        ZA1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[171][0])
        ZA1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[172][0])
        ZA1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[173][0])
        ZA1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[174][0])
        ZA2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[175][0])
        ZA2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[176][0])
        ZA2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[177][0])
        ZA2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[178][0])
        ZA2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[179][0])
        ZA2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[180][0])
        ZA2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[181][0])
        ZA12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[182][0])
        ZA12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[183][0])
        ZA12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[184][0])
        ZA12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[185][0])
        ZA12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[186][0])
        ZA12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[187][0])
        ZA12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[188][0])
        ZM1_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[189][0])
        ZM1_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[190][0])
        ZM1_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[191][0])
        ZM1_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[192][0])
        ZM1_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[193][0])
        ZM1_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[194][0])
        ZM1_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[195][0])
        ZM2_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[196][0])
        ZM2_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[197][0])
        ZM2_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[198][0])
        ZM2_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[199][0])
        ZM2_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[200][0])
        ZM2_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[201][0])
        ZM2_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[202][0])
        ZM12_MSR_MCMC_OPT.append(MCMC_CDF_VLS_OPT[203][0])
        ZM12_MSR_MC_1_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[204][0])
        ZM12_MSR_MC_1_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[205][0])
        ZM12_MSR_MC_2_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[206][0])
        ZM12_MSR_MC_2_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[207][0])
        ZM12_MSR_MC_3_SG_L_MCMC_OPT.append(MCMC_CDF_VLS_OPT[208][0])
        ZM12_MSR_MC_3_SG_H_MCMC_OPT.append(MCMC_CDF_VLS_OPT[209][0])

        
        Z1_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[0][0])
        Z2_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[2][0])
        Z1_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[4][0])
        Z2_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[6][0])
        V1_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[8][0])
        V2_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[10][0])
        V1_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[12][0])
        V2_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[14][0])
        S1_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[16][0])
        S2_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[18][0])
        S1_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[20][0])
        S2_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[22][0])
        L1_AVG_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[24][0])
        L1_MED_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[26][0])
        L1_AVG_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[28][0])
        L1_MED_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[30][0])
        L2_AVG_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[32][0])
        L2_MED_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[34][0])
        L2_AVG_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[36][0])
        L2_MED_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[38][0])
        S12_AVG_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[40][0])
        S12_MED_E1_MC_OPT.append(MCMC_CDF_PLT_OPT[42][0])
        L12_AVG_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[44][0])
        L12_MED_E1_MC_1_OPT.append(MCMC_CDF_PLT_OPT[46][0])
        L12_AVG_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[48][0])
        L12_MED_E1_MC_2_OPT.append(MCMC_CDF_PLT_OPT[50][0])

        Z1_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[1][0])
        Z2_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[3][0])
        Z1_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[5][0])
        Z2_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[7][0])
        V1_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[9][0])
        V2_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[11][0])
        V1_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[13][0])
        V2_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[15][0])
        S1_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[17][0])
        S2_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[19][0])
        S1_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[21][0])
        S2_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[23][0])
        L1_AVG_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[25][0])
        L1_MED_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[27][0])
        L1_AVG_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[29][0])
        L1_MED_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[31][0])
        L2_AVG_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[33][0])
        L2_MED_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[35][0])
        L2_AVG_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[37][0])
        L2_MED_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[39][0])
        S12_AVG_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[41][0])
        S12_MED_E2_MC_OPT.append(MCMC_CDF_PLT_OPT[43][0])
        L12_AVG_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[45][0])
        L12_MED_E2_MC_1_OPT.append(MCMC_CDF_PLT_OPT[47][0])
        L12_AVG_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[49][0])
        L12_MED_E2_MC_2_OPT.append(MCMC_CDF_PLT_OPT[51][0])

        #####################TABLE-1#####################
        mct1 = aptbl.Table()

        mct1['SMPL'] = SMPL

        mct1['FL1_AVG'] = FA1_MSR_MCMC_OPT
        mct1['FL1_AVG_MC_1SGL'] = FA1_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL1_AVG_MC_1SGH'] = FA1_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL1_AVG_MC_2SGL'] = FA1_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL1_AVG_MC_2SGH'] = FA1_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL1_AVG_MC_3SGL'] = FA1_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL1_AVG_MC_3SGH'] = FA1_MSR_MC_3_SG_H_MCMC_OPT

        mct1['FL2_AVG'] = FA2_MSR_MCMC_OPT
        mct1['FL2_AVG_MC_1SGL'] = FA2_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL2_AVG_MC_1SGH'] = FA2_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL2_AVG_MC_2SGL'] = FA2_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL2_AVG_MC_2SGH'] = FA2_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL2_AVG_MC_3SGL'] = FA2_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL2_AVG_MC_3SGH'] = FA2_MSR_MC_3_SG_H_MCMC_OPT

        mct1['FL12_AVG'] = FA12_MSR_MCMC_OPT
        mct1['FL12_AVG_MC_1SGL'] = FA12_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL12_AVG_MC_1SGH'] = FA12_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL12_AVG_MC_2SGL'] = FA12_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL12_AVG_MC_2SGH'] = FA12_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL12_AVG_MC_3SGL'] = FA12_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL12_AVG_MC_3SGH'] = FA12_MSR_MC_3_SG_H_MCMC_OPT

        mct1['FL1_MED'] = FM1_MSR_MCMC_OPT
        mct1['FL1_MED_MC_1SGL'] = FM1_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_1SGH'] = FM1_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL1_MED_MC_2SGL'] = FM1_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_2SGH'] = FM1_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL1_MED_MC_3SGL'] = FM1_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_3SGH'] = FM1_MSR_MC_3_SG_H_MCMC_OPT

        mct1['FL1_MED'] = FM2_MSR_MCMC_OPT
        mct1['FL1_MED_MC_1SGL'] = FM2_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_1SGH'] = FM2_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL1_MED_MC_2SGL'] = FM2_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_2SGH'] = FM2_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL1_MED_MC_3SGL'] = FM2_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL1_MED_MC_3SGH'] = FM2_MSR_MC_3_SG_H_MCMC_OPT

        mct1['FL12_MED'] = FM12_MSR_MCMC_OPT
        mct1['FL12_MED_MC_1SGL'] = FM12_MSR_MC_1_SG_L_MCMC_OPT
        mct1['FL12_MED_MC_1SGH'] = FM12_MSR_MC_1_SG_H_MCMC_OPT
        mct1['FL12_MED_MC_2SGL'] = FM12_MSR_MC_2_SG_L_MCMC_OPT
        mct1['FL12_MED_MC_2SGH'] = FM12_MSR_MC_2_SG_H_MCMC_OPT
        mct1['FL12_MED_MC_3SGL'] = FM12_MSR_MC_3_SG_L_MCMC_OPT
        mct1['FL12_MED_MC_3SGH'] = FM12_MSR_MC_3_SG_H_MCMC_OPT

        TABLESTATNAME_111  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_112  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'
        TABLESTATNAME_121  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_122  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct1.write(TABLESTATNAME_111, format='ascii.fixed_width_two_line', overwrite = True)
        mct1.write(TABLESTATNAME_112, format=tbl_format_opt, overwrite = True)
        mct1.write(TABLESTATNAME_121, format='ascii.fixed_width_two_line', overwrite = True)
        mct1.write(TABLESTATNAME_122, format=tbl_format_opt, overwrite = True)
        #####################TABLE-1#####################

        ###################TABLE-1-PLT###################
        mct1_plt = aptbl.Table()
        mct1_plt['SMPL'] = SMPL

        mct1_plt['Z1_AVG_E1_MC_OPT'] = Z1_AVG_E1_MC_OPT
        mct1_plt['Z2_AVG_E1_MC_OPT'] = Z2_AVG_E1_MC_OPT
        mct1_plt['Z1_MED_E1_MC_OPT'] = Z1_MED_E1_MC_OPT
        mct1_plt['Z2_MED_E1_MC_OPT'] = Z2_MED_E1_MC_OPT
        mct1_plt['V1_AVG_E1_MC_OPT'] = V1_AVG_E1_MC_OPT
        mct1_plt['V2_AVG_E1_MC_OPT'] = V2_AVG_E1_MC_OPT
        mct1_plt['V1_MED_E1_MC_OPT'] = V1_MED_E1_MC_OPT
        mct1_plt['V2_MED_E1_MC_OPT'] = V2_MED_E1_MC_OPT
        mct1_plt['S1_AVG_E1_MC_OPT'] = S1_AVG_E1_MC_OPT
        mct1_plt['S2_AVG_E1_MC_OPT'] = S2_AVG_E1_MC_OPT
        mct1_plt['S1_MED_E1_MC_OPT'] = S1_MED_E1_MC_OPT
        mct1_plt['S2_MED_E1_MC_OPT'] = S2_MED_E1_MC_OPT
        mct1_plt['L1_AVG_E1_MC_1_OPT'] = L1_AVG_E1_MC_1_OPT
        mct1_plt['L1_MED_E1_MC_1_OPT'] = L1_MED_E1_MC_1_OPT
        mct1_plt['L1_AVG_E1_MC_2_OPT'] = L1_AVG_E1_MC_2_OPT
        mct1_plt['L1_MED_E1_MC_2_OPT'] = L1_MED_E1_MC_2_OPT
        mct1_plt['L2_AVG_E1_MC_1_OPT'] = L2_AVG_E1_MC_1_OPT
        mct1_plt['L2_MED_E1_MC_1_OPT'] = L2_MED_E1_MC_1_OPT
        mct1_plt['L2_AVG_E1_MC_2_OPT'] = L2_AVG_E1_MC_2_OPT
        mct1_plt['L2_MED_E1_MC_2_OPT'] = L2_MED_E1_MC_2_OPT
        mct1_plt['S12_AVG_E1_MC_OPT'] = S12_AVG_E1_MC_OPT
        mct1_plt['S12_MED_E1_MC_OPT'] = S12_MED_E1_MC_OPT
        mct1_plt['L12_AVG_E1_MC_1_OPT'] = L12_AVG_E1_MC_1_OPT
        mct1_plt['L12_MED_E1_MC_1_OPT'] = L12_MED_E1_MC_1_OPT
        mct1_plt['L12_AVG_E1_MC_2_OPT'] = L12_AVG_E1_MC_2_OPT
        mct1_plt['L12_MED_E1_MC_2_OPT'] = L12_MED_E1_MC_2_OPT

        mct1_plt['Z1_AVG_E2_MC_OPT'] = Z1_AVG_E2_MC_OPT
        mct1_plt['Z2_AVG_E2_MC_OPT'] = Z2_AVG_E2_MC_OPT
        mct1_plt['Z1_MED_E2_MC_OPT'] = Z1_MED_E2_MC_OPT
        mct1_plt['Z2_MED_E2_MC_OPT'] = Z2_MED_E2_MC_OPT
        mct1_plt['V1_AVG_E2_MC_OPT'] = V1_AVG_E2_MC_OPT
        mct1_plt['V2_AVG_E2_MC_OPT'] = V2_AVG_E2_MC_OPT
        mct1_plt['V1_MED_E2_MC_OPT'] = V1_MED_E2_MC_OPT
        mct1_plt['V2_MED_E2_MC_OPT'] = V2_MED_E2_MC_OPT
        mct1_plt['S1_AVG_E2_MC_OPT'] = S1_AVG_E2_MC_OPT
        mct1_plt['S2_AVG_E2_MC_OPT'] = S2_AVG_E2_MC_OPT
        mct1_plt['S1_MED_E2_MC_OPT'] = S1_MED_E2_MC_OPT
        mct1_plt['S2_MED_E2_MC_OPT'] = S2_MED_E2_MC_OPT
        mct1_plt['L1_AVG_E2_MC_1_OPT'] = L1_AVG_E2_MC_1_OPT
        mct1_plt['L1_MED_E2_MC_1_OPT'] = L1_MED_E2_MC_1_OPT
        mct1_plt['L1_AVG_E2_MC_2_OPT'] = L1_AVG_E2_MC_2_OPT
        mct1_plt['L1_MED_E2_MC_2_OPT'] = L1_MED_E2_MC_2_OPT
        mct1_plt['L2_AVG_E2_MC_1_OPT'] = L2_AVG_E2_MC_1_OPT
        mct1_plt['L2_MED_E2_MC_1_OPT'] = L2_MED_E2_MC_1_OPT
        mct1_plt['L2_AVG_E2_MC_2_OPT'] = L2_AVG_E2_MC_2_OPT
        mct1_plt['L2_MED_E2_MC_2_OPT'] = L2_MED_E2_MC_2_OPT
        mct1_plt['S12_AVG_E2_MC_OPT'] = S12_AVG_E2_MC_OPT
        mct1_plt['S12_MED_E2_MC_OPT'] = S12_MED_E2_MC_OPT
        mct1_plt['L12_AVG_E2_MC_1_OPT'] = L12_AVG_E2_MC_1_OPT
        mct1_plt['L12_MED_E2_MC_1_OPT'] = L12_MED_E2_MC_1_OPT
        mct1_plt['L12_AVG_E2_MC_2_OPT'] = L12_AVG_E2_MC_2_OPT
        mct1_plt['L12_MED_E2_MC_2_OPT'] = L12_MED_E2_MC_2_OPT

        TABLESTATNAME_11_PLT  = plt_dir_tbl + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt-PLT.dat'
        TABLESTATNAME_12_PLT  = plt_dir_tbl + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt-PLT.csv'
        TABLESTATNAME_21_PLT  = plt_dir_tbl + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt-PLT.dat'
        TABLESTATNAME_22_PLT  = plt_dir_tbl + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) +'-FLX-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt-PLT.csv'

        mct1_plt.write(TABLESTATNAME_11_PLT, format='ascii.fixed_width_two_line', overwrite = True)
        mct1_plt.write(TABLESTATNAME_12_PLT, format=tbl_format_opt, overwrite = True)
        mct1_plt.write(TABLESTATNAME_21_PLT, format='ascii.fixed_width_two_line', overwrite = True)
        mct1_plt.write(TABLESTATNAME_22_PLT, format=tbl_format_opt, overwrite = True)

        ###################TABLE-1-PLT###################       

        #####################TABLE-2#####################
        mct2 = aptbl.Table()
        mct2['SMPL'] = SMPL

        mct2['LM1_AVG'] = LA1_MSR_MCMC_OPT
        mct2['LM1_AVG_MC_1SGL'] = LA1_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM1_AVG_MC_1SGH'] = LA1_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM1_AVG_MC_2SGL'] = LA1_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM1_AVG_MC_2SGH'] = LA1_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM1_AVG_MC_3SGL'] = LA1_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM1_AVG_MC_3SGH'] = LA1_MSR_MC_3_SG_H_MCMC_OPT

        mct2['LM2_AVG'] = LA2_MSR_MCMC_OPT
        mct2['LM2_AVG_MC_1SGL'] = LA2_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM2_AVG_MC_1SGH'] = LA2_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM2_AVG_MC_2SGL'] = LA2_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM2_AVG_MC_2SGH'] = LA2_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM2_AVG_MC_3SGL'] = LA2_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM2_AVG_MC_3SGH'] = LA2_MSR_MC_3_SG_H_MCMC_OPT

        mct2['LM12_AVG'] = LA12_MSR_MCMC_OPT
        mct2['LM12_AVG_MC_1SGL'] = LA12_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM12_AVG_MC_1SGH'] = LA12_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM12_AVG_MC_2SGL'] = LA12_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM12_AVG_MC_2SGH'] = LA12_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM12_AVG_MC_3SGL'] = LA12_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM12_AVG_MC_3SGH'] = LA12_MSR_MC_3_SG_H_MCMC_OPT

        mct2['LM1_MED'] = LM1_MSR_MCMC_OPT
        mct2['LM1_MED_MC_1SGL'] = LM1_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM1_MED_MC_1SGH'] = LM1_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM1_MED_MC_2SGL'] = LM1_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM1_MED_MC_2SGH'] = LM1_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM1_MED_MC_3SGL'] = LM1_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM1_MED_MC_3SGH'] = LM1_MSR_MC_3_SG_H_MCMC_OPT

        mct2['LM2_MED'] = LM2_MSR_MCMC_OPT
        mct2['LM2_MED_MC_1SGL'] = LM2_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM2_MED_MC_1SGH'] = LM2_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM2_MED_MC_2SGL'] = LM2_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM2_MED_MC_2SGH'] = LM2_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM2_MED_MC_3SGL'] = LM2_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM2_MED_MC_3SGH'] = LM2_MSR_MC_3_SG_H_MCMC_OPT

        mct2['LM12_MED'] = LM12_MSR_MCMC_OPT
        mct2['LM12_MED_MC_1SGL'] = LM12_MSR_MC_1_SG_L_MCMC_OPT
        mct2['LM12_MED_MC_1SGH'] = LM12_MSR_MC_1_SG_H_MCMC_OPT
        mct2['LM12_MED_MC_2SGL'] = LM12_MSR_MC_2_SG_L_MCMC_OPT
        mct2['LM12_MED_MC_2SGH'] = LM12_MSR_MC_2_SG_H_MCMC_OPT
        mct2['LM12_MED_MC_3SGL'] = LM12_MSR_MC_3_SG_L_MCMC_OPT
        mct2['LM12_MED_MC_3SGH'] = LM12_MSR_MC_3_SG_H_MCMC_OPT

        TABLESTATNAME_211  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_212  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'
        TABLESTATNAME_221  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_222  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct2.write(TABLESTATNAME_211, format='ascii.fixed_width_two_line', overwrite = True)
        mct2.write(TABLESTATNAME_212, format=tbl_format_opt, overwrite = True)
        mct2.write(TABLESTATNAME_221, format='ascii.fixed_width_two_line', overwrite = True)
        mct2.write(TABLESTATNAME_222, format=tbl_format_opt, overwrite = True)
        #####################TABLE-2#####################

        #####################TABLE-3#####################
        mct3 = aptbl.Table()
        mct3['SMPL'] = SMPL

        mct3['LLM1_AVG'] = LLA1_MSR_MCMC_OPT
        mct3['LLM1_AVG_MC_1SGL'] = LLA1_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM1_AVG_MC_1SGH'] = LLA1_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM1_AVG_MC_2SGL'] = LLA1_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM1_AVG_MC_2SGH'] = LLA1_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM1_AVG_MC_3SGL'] = LLA1_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM1_AVG_MC_3SGH'] = LLA1_MSR_MC_3_SG_H_MCMC_OPT

        mct3['LLM2_AVG'] = LLA2_MSR_MCMC_OPT
        mct3['LLM2_AVG_MC_1SGL'] = LLA2_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM2_AVG_MC_1SGH'] = LLA2_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM2_AVG_MC_2SGL'] = LLA2_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM2_AVG_MC_2SGH'] = LLA2_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM2_AVG_MC_3SGL'] = LLA2_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM2_AVG_MC_3SGH'] = LLA2_MSR_MC_3_SG_H_MCMC_OPT

        mct3['LLM12_AVG'] = LLA12_MSR_MCMC_OPT
        mct3['LLM12_AVG_MC_1SGL'] = LLA12_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM12_AVG_MC_1SGH'] = LLA12_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM12_AVG_MC_2SGL'] = LLA12_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM12_AVG_MC_2SGH'] = LLA12_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM12_AVG_MC_3SGL'] = LLA12_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM12_AVG_MC_3SGH'] = LLA12_MSR_MC_3_SG_H_MCMC_OPT

        mct3['LLM1_MED'] = LLM1_MSR_MCMC_OPT
        mct3['LLM1_MED_MC_1SGL'] = LLM1_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM1_MED_MC_1SGH'] = LLM1_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM1_MED_MC_2SGL'] = LLM1_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM1_MED_MC_2SGH'] = LLM1_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM1_MED_MC_3SGL'] = LLM1_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM1_MED_MC_3SGH'] = LLM1_MSR_MC_3_SG_H_MCMC_OPT

        mct3['LLM2_MED'] = LLM2_MSR_MCMC_OPT
        mct3['LLM2_MED_MC_1SGL'] = LLM2_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM2_MED_MC_1SGH'] = LLM2_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM2_MED_MC_2SGL'] = LLM2_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM2_MED_MC_2SGH'] = LLM2_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM2_MED_MC_3SGL'] = LLM2_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM2_MED_MC_3SGH'] = LLM2_MSR_MC_3_SG_H_MCMC_OPT

        mct3['LLM12_MED'] = LLM12_MSR_MCMC_OPT
        mct3['LLM12_MED_MC_1SGL'] = LLM12_MSR_MC_1_SG_L_MCMC_OPT
        mct3['LLM12_MED_MC_1SGH'] = LLM12_MSR_MC_1_SG_H_MCMC_OPT
        mct3['LLM12_MED_MC_2SGL'] = LLM12_MSR_MC_2_SG_L_MCMC_OPT
        mct3['LLM12_MED_MC_2SGH'] = LLM12_MSR_MC_2_SG_H_MCMC_OPT
        mct3['LLM12_MED_MC_3SGL'] = LLM12_MSR_MC_3_SG_L_MCMC_OPT
        mct3['LLM12_MED_MC_3SGH'] = LLM12_MSR_MC_3_SG_H_MCMC_OPT

        TABLESTATNAME_311  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_312  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'
        TABLESTATNAME_321  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_322  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-LUM-LOG-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct3.write(TABLESTATNAME_311, format='ascii.fixed_width_two_line', overwrite = True)
        mct3.write(TABLESTATNAME_312, format=tbl_format_opt, overwrite = True)
        mct3.write(TABLESTATNAME_321, format='ascii.fixed_width_two_line', overwrite = True)
        mct3.write(TABLESTATNAME_322, format=tbl_format_opt, overwrite = True)
        #####################TABLE-3#####################

        #####################TABLE-4#####################
        mct4 = aptbl.Table()
        mct4['SMPL'] = SMPL

        mct4['VA1_AVG'] = VA1_MSR_MCMC_OPT
        mct4['VA1_AVG_MC_1SGL'] = VA1_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA1_AVG_MC_1SGH'] = VA1_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA1_AVG_MC_2SGL'] = VA1_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA1_AVG_MC_2SGH'] = VA1_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA1_AVG_MC_3SGL'] = VA1_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA1_AVG_MC_3SGH'] = VA1_MSR_MC_3_SG_H_MCMC_OPT

        mct4['VA2_AVG'] = VA2_MSR_MCMC_OPT
        mct4['VA2_AVG_MC_1SGL'] = VA2_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA2_AVG_MC_1SGH'] = VA2_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA2_AVG_MC_2SGL'] = VA2_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA2_AVG_MC_2SGH'] = VA2_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA2_AVG_MC_3SGL'] = VA2_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA2_AVG_MC_3SGH'] = VA2_MSR_MC_3_SG_H_MCMC_OPT

        mct4['VA12_AVG'] = VA12_MSR_MCMC_OPT
        mct4['VA12_AVG_MC_1SGL'] = VA12_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA12_AVG_MC_1SGH'] = VA12_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA12_AVG_MC_2SGL'] = VA12_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA12_AVG_MC_2SGH'] = VA12_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA12_AVG_MC_3SGL'] = VA12_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA12_AVG_MC_3SGH'] = VA12_MSR_MC_3_SG_H_MCMC_OPT

        mct4['VM1_MED'] = VM1_MSR_MCMC_OPT
        mct4['VA1_MED_MC_1SGL'] = VM1_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA1_MED_MC_1SGH'] = VM1_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA1_MED_MC_2SGL'] = VM1_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA1_MED_MC_2SGH'] = VM1_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA1_MED_MC_3SGL'] = VM1_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA1_MED_MC_3SGH'] = VM1_MSR_MC_3_SG_H_MCMC_OPT

        mct4['VA2_MED'] = VM2_MSR_MCMC_OPT
        mct4['VA2_MED_MC_1SGL'] = VM2_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA2_MED_MC_1SGH'] = VM2_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA2_MED_MC_2SGL'] = VM2_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA2_MED_MC_2SGH'] = VM2_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA2_MED_MC_3SGL'] = VM2_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA2_MED_MC_3SGH'] = VM2_MSR_MC_3_SG_H_MCMC_OPT

        mct4['VA12_MED'] = VM12_MSR_MCMC_OPT
        mct4['VA12_MED_MC_1SGL'] = VM12_MSR_MC_1_SG_L_MCMC_OPT
        mct4['VA12_MED_MC_1SGH'] = VM12_MSR_MC_1_SG_H_MCMC_OPT
        mct4['VA12_MED_MC_2SGL'] = VM12_MSR_MC_2_SG_L_MCMC_OPT
        mct4['VA12_MED_MC_2SGH'] = VM12_MSR_MC_2_SG_H_MCMC_OPT
        mct4['VA12_MED_MC_3SGL'] = VM12_MSR_MC_3_SG_L_MCMC_OPT
        mct4['VA12_MED_MC_3SGH'] = VM12_MSR_MC_3_SG_H_MCMC_OPT

        TABLESTATNAME_411  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_412  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'
        TABLESTATNAME_421  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_422  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-'+str(sbsms_mcmc)+'-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct4.write(TABLESTATNAME_411, format='ascii.fixed_width_two_line', overwrite = True)
        mct4.write(TABLESTATNAME_412, format=tbl_format_opt, overwrite = True)
        mct4.write(TABLESTATNAME_421, format='ascii.fixed_width_two_line', overwrite = True)
        mct4.write(TABLESTATNAME_422, format=tbl_format_opt, overwrite = True)
        #####################TABLE-4#####################

        #####################TABLE-5#####################
        mct5 = aptbl.Table()
        mct5['SMPL'] = SMPL

        mct5['z1_AVG'] = ZA1_MSR_MCMC_OPT
        mct5['z1_MSR_MC_1SGL'] = ZA1_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z1_MSR_MC_1SGH'] = ZA1_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z1_MSR_MC_2SGL'] = ZA1_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z1_MSR_MC_2SGH'] = ZA1_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z1_MSR_MC_3SGL'] = ZA1_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z1_MSR_MC_3SGH'] = ZA1_MSR_MC_3_SG_H_MCMC_OPT

        mct5['z2_AVG'] = ZA2_MSR_MCMC_OPT
        mct5['z2_AVG_MC_1SGL'] = ZA2_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z2_AVG_MC_1SGH'] = ZA2_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z2_AVG_MC_2SGL'] = ZA2_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z2_AVG_MC_2SGH'] = ZA2_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z2_AVG_MC_3SGL'] = ZA2_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z2_AVG_MC_3SGH'] = ZA2_MSR_MC_3_SG_H_MCMC_OPT

        mct5['z12_AVG'] = ZA12_MSR_MCMC_OPT
        mct5['z12_AVG_MC_1SGL'] = ZA12_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z12_AVG_MC_1SGH'] = ZA12_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z12_AVG_MC_2SGL'] = ZA12_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z12_AVG_MC_2SGH'] = ZA12_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z12_AVG_MC_3SGL'] = ZA12_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z12_AVG_MC_3SGH'] = ZA12_MSR_MC_3_SG_H_MCMC_OPT

        mct5['z1_MED'] = ZM1_MSR_MCMC_OPT
        mct5['z1_MED_MC_1SGL'] = ZM1_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z1_MED_MC_1SGH'] = ZM1_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z1_MED_MC_2SGL'] = ZM1_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z1_MED_MC_2SGH'] = ZM1_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z1_MED_MC_3SGL'] = ZM1_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z1_MED_MC_3SGH'] = ZM1_MSR_MC_3_SG_H_MCMC_OPT

        mct5['z2_MED'] = ZM2_MSR_MCMC_OPT
        mct5['z2_MED_MC_1SGL'] = ZM2_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z2_MED_MC_1SGH'] = ZM2_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z2_MED_MC_2SGL'] = ZM2_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z2_MED_MC_2SGH'] = ZM2_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z2_MED_MC_3SGL'] = ZM2_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z2_MED_MC_3SGH'] = ZM2_MSR_MC_3_SG_H_MCMC_OPT

        mct5['z12_MED'] = ZM12_MSR_MCMC_OPT
        mct5['z12_MED_MC_1SGL'] = ZM12_MSR_MC_1_SG_L_MCMC_OPT
        mct5['z12_MED_MC_1SGH'] = ZM12_MSR_MC_1_SG_H_MCMC_OPT
        mct5['z12_MED_MC_2SGL'] = ZM12_MSR_MC_2_SG_L_MCMC_OPT
        mct5['z12_MED_MC_2SGH'] = ZM12_MSR_MC_2_SG_H_MCMC_OPT
        mct5['z12_MED_MC_3SGL'] = ZM12_MSR_MC_3_SG_L_MCMC_OPT
        mct5['z12_MED_MC_3SGH'] = ZM12_MSR_MC_3_SG_H_MCMC_OPT

        TABLESTATNAME_511  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_512  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'
        TABLESTATNAME_521  = mcm_dir_plt_1 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.dat'
        TABLESTATNAME_522  = mcm_dir_plt_2 + '/CII_HATLAS-' + line1 + '-' + line2+'-' +str(sbsms_mcmc) + '-MS-' +str(method) + '-MC-'+str(iterations_mc)+'-' +str(method) + '-Z-'+str(sbsmn_mcmc[0])+'-'+str(sbsmn_mcmc[-1])+'-stk-'+str(spc_wdt_dir)+'kms-crc-'+str(mask_radi_as_ms)+'as_msk_ms-stt.csv'

        mct5.write(TABLESTATNAME_511, format='ascii.fixed_width_two_line', overwrite = True)
        mct5.write(TABLESTATNAME_512, format=tbl_format_opt, overwrite = True)
        mct5.write(TABLESTATNAME_521, format='ascii.fixed_width_two_line', overwrite = True)
        mct5.write(TABLESTATNAME_522, format=tbl_format_opt, overwrite = True)
        #####################TABLE-5#####################

    print
    print (colored('MC Statistics table: '+TABLESTATNAME_1_1_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_1_1_2,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_1_2_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_1_2_2,'magenta'))
    print
    print (colored('MC Statistics table: '+TABLESTATNAME_2_1_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_2_1_2,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_2_2_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_2_2_2,'magenta'))
    print
    print (colored('MC Statistics table: '+TABLESTATNAME_3_1_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_3_1_2,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_3_2_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_3_2_2,'magenta'))

    print
    print (colored('MC Statistics table: '+TABLESTATNAME_4_1_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_4_1_2,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_4_2_1,'magenta'))
    print (colored('MC Statistics table: '+TABLESTATNAME_4_2_2,'magenta'))

    print
    print (colored('MC 1 Statistics table: '+TABLESTATNAME_111,'green'))
    print (colored('MC 1 Statistics table: '+TABLESTATNAME_112,'green'))
    print (colored('MC 1 Statistics table: '+TABLESTATNAME_121,'green'))
    print (colored('MC 1 Statistics table: '+TABLESTATNAME_122,'green'))
    print
    print (colored('MC 2 Statistics table: '+TABLESTATNAME_211,'green'))
    print (colored('MC 2 Statistics table: '+TABLESTATNAME_212,'green'))
    print (colored('MC 2 Statistics table: '+TABLESTATNAME_221,'green'))
    print (colored('MC 2 Statistics table: '+TABLESTATNAME_222,'green'))
    print
    print (colored('MC 3 Statistics table: '+TABLESTATNAME_311,'green'))
    print (colored('MC 3 Statistics table: '+TABLESTATNAME_312,'green'))
    print (colored('MC 3 Statistics table: '+TABLESTATNAME_321,'green'))
    print (colored('MC 3 Statistics table: '+TABLESTATNAME_322,'green'))
    print
    print (colored('MC 4 Statistics table: '+TABLESTATNAME_411,'green'))
    print (colored('MC 4 Statistics table: '+TABLESTATNAME_412,'green'))
    print (colored('MC 4 Statistics table: '+TABLESTATNAME_421,'green'))
    print (colored('MC 4 Statistics table: '+TABLESTATNAME_422,'green'))
    print
    print (colored('MC 5 Statistics table: '+TABLESTATNAME_511,'green'))
    print (colored('MC 5 Statistics table: '+TABLESTATNAME_512,'green'))
    print (colored('MC 5 Statistics table: '+TABLESTATNAME_521,'green'))
    print (colored('MC 5 Statistics table: '+TABLESTATNAME_522,'green'))
    print

    print
    print (colored('MC 1 Statistics table (Plots): ' + TABLESTATNAME_11_PLT,'green'))
    print (colored('MC 1 Statistics table (Plots): ' + TABLESTATNAME_12_PLT,'green'))
    print (colored('MC 1 Statistics table (Plots): ' + TABLESTATNAME_21_PLT,'green'))
    print (colored('MC 1 Statistics table (Plots): ' + TABLESTATNAME_22_PLT,'green'))
    print
####Fnc_Stk_Plt###