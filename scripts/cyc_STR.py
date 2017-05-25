from netCDF4 import Dataset
import numpy as np
import datetime as datetime
import calendar
import matplotlib.pyplot as plt
from scipy import stats
# ************************************************
#  Selection
# ************************************************
month_abbr = ["","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"]

# time
year = [1979,1984]
yrs = np.arange(year[0],year[1],1)
nyrs = np.size(yrs)
x = np.arange(0,nyrs-1,1)

mm = [1,12]  #[ 9,9]
dd = [1,31]
hr = [0,18]

# ssn   / use "YYY" for whole year
ssn="YYY"

# nf
ntrk1 = 0# 5518 ;0
ntrk2 = 20000# 5518 ;20000

hs = "SH"
# reg

lon = [ -90,361]
lat = [ -40,-20]
reg = "20_40S"

# lon = [ 0,150]
# lat = [ -40,-20]
# reg = "0_150E.20_40S"

maxnf = 200 # max # of fronts per timestep
maxnp = 100 # max # of frontal points

# ************************************************
# read in netCDF files
# ************************************************
# 1. slp
# ************************************************

# for yr in yrs:
#   slpyr = yr
#   fin = "/Users/irudeva/work/DATA/ERAint/Mslp_highres/erain.mslp.avm.%d.nc"%slpyr
#   print fin
#   nc = Dataset(fin, 'r')
#
#   # varnam=['longitude','latitude','time','msl']
#   # v=0
#   # for var in varnam:
#   #   if nc.variables[varnam[v]].name != var:
#   #       print "Variables don't agree", var, nc.variables[varnam[v]].name, v
#   #       exit()
#   #   v += 1
#
#   # nc_attrs, nc_dims, nc_vars = ncdump(nc)
#
#   # print nc_attrs, nc_dims, nc_vars
#   print nc.variables.keys()
#   if nc.variables.keys()[0]=="lon":
#       lonslp = nc.variables['lon'][:]
#       latslp = nc.variables['lat'][:]
#   elif nc.variables.keys()[0]=="longitude":
#       lonslp = nc.variables['longitude'][:]
#       latslp = nc.variables['lattitude'][:]
#   else:
#      print"Check lon/lat names in the netcdf file, %d "%yr
#   time = nc.variables['time'][:]
#   mslp = nc.variables['msl'][:]/100.
#
#   dt_nc = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
#        for t in time]
#   print dt_nc[1].year
#
#
#   latslpSH = latslp[np.logical_and(latslp<-20,latslp>-65)]
#   mslpSH   = mslp[:,np.logical_and(latslp<-20,latslp>-65),:]
#
#   # zonal average
#   slp_z = np.mean(mslpSH[:,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=2)
#
#   if yr==year[0]:
#       STRslp = np.zeros((2,12,nyrs),dtype=np.float)
#       STRlat = np.zeros_like(STRslp)
#
#   for im in np.arange(slp_z[:,0].size):
#     #   print dt_nc[im].month
#       STRslp[0,im,yr-year[0]] = np.amax(slp_z[im,:])
#       STRlat[0,im,yr-year[0]] = latslpSH[np.argmax(slp_z[im,:])]
#     #   print STRlat[im,yr-year[0]], STRslp[im,yr-year[0]]
#
#   nc.close()
#

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Plotting
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plt.close('all')
#
# f, ax = plt.subplots(6, 2,  sharex='col', sharey='row')
# plt.suptitle('STR intensity, SH',fontsize=14)
# for ic in np.arange(2):
#     for ir in np.arange(6):
#         im = ir+ic*6
#         a = ax[ir, ic]
#         a2 = a.twinx()
#
#         a.plot(yrs,STRslp[im,:],color='b',lw=1,label='STR intensity')
#         a2.plot(yrs,STRlat[im,:],color='g',lw=1,label='STR location')
#         # ax[ir, ic].set_title('STR intensity, %s, SH'%month_abbr[ir+ic*6+1])
#         a.axis((year[0]-1, year[1], 1010, 1025))
#         a2.axis((year[0]-1, year[1], -40, -25))
#
#         if im == 3:
#             ax[ir, ic].set_ylabel('MSL Pressure (Pa)',color='b')
#         if im == 9:
#             a2.set_ylabel('latitude',color='g')
#
#         #  trend
#         slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STRslp[im,:])
#         trend = intercept + (slope * yrs)
#         a.plot(yrs,trend,color='b',lw=2)
#
#         slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STRlat[im,:])
#         trend = intercept + (slope * yrs)
#         a2.plot(yrs,trend,color='g',lw=2)
#
#         a.set_title(' %s, SH'%month_abbr[ir+ic*6+1])
#
#         # ax2[ir, ic] = ax[ir, ic].twinx()
#         # ax2.plot(t, s2, 'r.')
#         # ax2.set_ylabel('sin', color='r')
#
#         if im==12:
#             ax[ir, ic].legend(loc='lower right')
#
# plt.setp([a.set_xlabel('Year') for a in ax[5, :]])
# a=ax[-1, -1]
# # ax[0, 0].xlabel('Year')
# # ax[0, 0].ylabel('MSL Pressure (Pa)')
# # ax[0, 0].axis((1979, 2014, 1012, 1030))
#
#
# f.subplots_adjust(hspace=0.3)
#
# plt.show()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Fronts
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for yr in yrs:
  fin = "/Users/Irina/work/Projects/trkUOM/ERAint/trkgrid/trk.%d.nc"%yr
  print "cyclones from %s"%fin
  ncf = Dataset(fin, 'r')

  print ncf.variables.keys()
  trklon = ncf.variables['trklon'][:]
  trklat = ncf.variables['trklat'][:]
  trktime = ncf.variables['trktime'][:]
  trkslp = ncf.variables['trkslp'][:]
  trklen = ncf.variables['trklen'][:]

  # fdt = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
  #      for t in time]

  #
  # if any(x > maxnf) :
  #   print "ERROR: nf > maxnf"
  #   quit()
  #
  # if any(x > maxnp ) :
  #   print "ERROR:  np > maxnp"
  #   quit()
  #
  # ************************************************
  # masking
  # ************************************************
  # create mask
  trkmask = np.zeros([nyrs,trklen.size*2,200],dtype=np.int)
  trkslp_mask = np.zeros_like(trkmask)

  for ntrk in np.arange(trklen.size):
      if ntrk+1>=ntrk1 and ntrk+1<=ntrk2:

             nit = trklen[ntrk]
             trkdate = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
                  for t in trktime[ntrk,:nit]]
             for n in np.arange(nit):
                 trkyr = trkdate[n].year
                 trkmon = trkdate[n].month
                 trkday = trkdate[n].day
                 trkhr = trkdate[n].hour

                 if trkmon in [12,1,2] :
                     trkssn = "DJF"
                 if ssn == "DJF" and trkmon ==12:
                     trkyr = trkyr+1
                 if trkmon in [3,4,5]:
                     trkssn = "MAM"
                 if trkmon in [6,7,8]:
                     trkssn = "JJA"
                 if trkmon in [9,10,11]:
                     trkssn = "SOP"


                 if trkyr>=year[0] and trkyr<=year[1]:
                     if trkmon>=mm[0] and trkmon<=mm[1]:
                         if trkday>=dd[0] and trkday<=dd[1]:
                             if trkhr>=hr[0] and trkhr<=hr[1]:

                                 if ssn == trkssn:

                                   if  trklat[ntrk,n]>=lat[0] and trklat[ntrk,n]<=lat[1]:
                                       if trklon[ntrk,n]>=lon[0] and trklon[ntrk,n]<=lon[1]:
                                           trkmask[yr-yrar[0],ntrk,n] = trkmon
                                           trkslp_mask[yr-yrar[0],ntrk,n] = trkslp[ntrk,n]



# ************************************************
# cyc stats
# ************************************************

crit = np.zeros(12,dtype = np.float)
ncyc_my    = np.zeros([ 2,12,nyrs],dtype=np.float)
ncyc_my3tr = np.zeros_like(ncyc_my)
# fr_northlat_my = np.zeros_like(ncyc_my)
# fr_northlat_my3tr = np.zeros_like(ncyc_my)

for im in range( 1,13):

  trkslptmp = trkslp_mask[np.where(trkmask == im)]

  crit[im-1] =  np.percentile(trkslptmp,67)
  print  month_abbr[im], crit[im-1]

  # hist, bin_edges = np.histogram(trlslptmp,range=(0.1,25))
  # print "pdf: ",hist
  # print bin_edges


for yr in yrs:
    for im in range( 1,13):
        ncyc_my[0,im-1,yr-year[0]]=np.sum(trkmask[yr-year[0],:,:]==im)

        # tmp0 = fr_northlat[yr-year[0],:,:][np.where(trkmask[yr-year[0],:,:] == im)]
        # fr_northlat_my[0,im-1,yr-year[0]]=np.mean(tmp0)

        tmp  = trlslp[yr-year[0],:,:][np.where(trkmask[yr-year[0],:,:] == im)]
        ncyc_my3tr[0,im-1,yr-year[0]]=np.count_nonzero(np.where(tmp>=crit[im-1], 1,0))
        # fr_northlat_my3tr[0,im-1,yr-year[0]]=np.mean(tmp0[np.where(tmp>=crit[im-1])])
        # # ncyc_my3tr[im-1,yr-year[0]]=np.sum(trkmask[yr-year[0],:,:]==im and trlslp[yr-year[0],:,:]>=crit[im-1])
        # # print yr, im, month_abbr[im], fr_northlat_my[0,im-1,yr-year[0]],fr_northlat_my3tr[0,im-1,yr-year[0]]
        # # print np.count_nonzero(np.where(trkmask[yr-year[0],:,:]== im, 1,0))

# ************************************************
# trends
# ************************************************

for im in range( 1,13):
    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,ncyc_my[0,im-1,:])
    ncyc_my[1,im-1,:] = intercept + (slope * yrs)

    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,ncyc_my3tr[0,im-1,:])
    ncyc_my3tr[1,im-1,:] = intercept + (slope * yrs)

    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,fr_northlat_my[0,im-1,:])
    fr_northlat_my[1,im-1,:] = intercept + (slope * yrs)

    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,fr_northlat_my3tr[0,im-1,:])
    fr_northlat_my3tr[1,im-1,:] = intercept + (slope * yrs)

    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STRlat[0,im-1,:])
    STRlat[1,im-1,:] = intercept + (slope * yrs)

    slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STRslp[0,im-1,:])
    STRslp[1,im-1,:] = intercept + (slope * yrs)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.close('all')
for ifig in np.arange(4,8):
    # fig = plt.figure(ifig)
    f, ax = plt.subplots(6, 2,  sharex='col', sharey='row')

    if ifig ==0:
        plt.suptitle('STR intensity vs number of fronts, %s'%reg,fontsize=14)
        fout = "../output/frN_STRint.%s.png"%reg
    if ifig ==1:
        plt.suptitle('STR intensity vs number of strong fronts, %s'%reg,fontsize=14)
        fout = "../output/frNstr_STRint.%s.png"%reg
    if ifig ==2:
        plt.suptitle('STR location vs number of fronts, %s'%reg,fontsize=14)
        fout = "../output/frN_STRloc.%s.png"%reg
    if ifig ==3:
        plt.suptitle('STR location vs number of strong fronts, %s'%reg,fontsize=14)
        fout = "../output/frNstr_STRloc.%s.png"%reg
    if ifig ==4:
        plt.suptitle('STR intensity vs north lat of fronts, %s'%reg,fontsize=14)
        fout = "../output/frNorthLat_STRint.%s.png"%reg
    if ifig ==5:
        plt.suptitle('STR intensity vs north lat of strong fronts, %s'%reg,fontsize=14)
        fout = "../output/frNorthLatstr_STRint.%s.png"%reg
    if ifig ==6:
        plt.suptitle('STR location vs north lat of fronts, %s'%reg,fontsize=14)
        fout = "../output/frNorthLat_STRloc.%s.png"%reg
    if ifig ==7:
        plt.suptitle('STR location vs north lat of strong fronts, %s'%reg,fontsize=14)
        fout = "../output/frNorthLatstr_STRloc.%s.png"%reg

    if ifig ==0 or ifig == 2:
        var1 = ncyc_my
        yax1 = [900, 1800]
        yax1label = 'number of fronts'
        var1legend = 'N of fronts'

    if ifig ==1 or ifig == 3:
        var1 = ncyc_my3tr
        yax1 = [300, 700]
        yax1label = 'number of fronts'
        var1legend = 'N of strong fronts (int > 67th perc)'

    if any ([4,6] == ifig ):
        var1 = fr_northlat_my
        yax1 = [-40, -25]
        yax1label = 'latitude'
        var1legend = 'front lat north'

    if any ([5,7] == ifig ):
        var1 = fr_northlat_my3tr
        yax1 = [-40, -25]
        yax1label = 'latitude'
        var1legend = 'strong front lat north'


    if any ([0,1,4,5] == ifig ):
        var2 = STRslp
        yax2 = [1010, 1025]
        yax2label = 'MSL Pressure (Pa)'
        var2legend = 'STR intendity'

    if any ([2,3,6,7] == ifig ):
        var2 = STRlat
        yax2 = [-40, -25]
        yax2label = 'latitude'
        var2legend = 'STR location'


    for ic in np.arange(2):
        for ir in np.arange(6):
            im = ir+ic*6
            a = ax[ir, ic]
            a2 = a.twinx()

            a.plot(yrs,var1[0,im,:],color='b',lw=1,label=var1legend)
            a2.plot(yrs,var2[0,im,:],color='g',lw=1,label=var2legend)
            # ax[ir, ic].set_title('STR intensity, %s, SH'%month_abbr[ir+ic*6+1])
            a.axis((year[0]-1, year[1], yax1[0], yax1[1]))
            a2.axis((year[0]-1, year[1], yax2[0], yax2[1]))

            if im == 3:
                ax[ir, ic].set_ylabel(yax1label,color='b')
            if im == 9:
                a2.set_ylabel(yax2label,color='g')

            #  trend
            a.plot(yrs,var1[1,im,:],color='b',lw=2)
            a2.plot(yrs,var2[1,im,:],color='g',lw=2)

            cc = np.corrcoef(var1[0,im,:],var2[0,im,:])
            cc_dt = np.corrcoef(var1[0,im,:]-var1[1,im,:],var2[0,im,:]-var2[1,im,:])
            a.set_title(' %s, %s, r = %.2f, r_dt = %.2f'%(month_abbr[ir+ic*6+1],reg,cc[1,0],cc_dt[1,0]))

            if im==11:
                ax[ir, ic].legend(loc='lower right')

    plt.setp([a.set_xlabel('Year') for a in ax[5, :]])
    # a=ax[-1, -1]
    # ax[0, 0].xlabel('Year')
    # ax[0, 0].ylabel('MSL Pressure (Pa)')
    # ax[0, 0].axis((1979, 2014, 1012, 1030))


    f.subplots_adjust(hspace=0.3)

    plt.show()
    f.savefig(fout)
