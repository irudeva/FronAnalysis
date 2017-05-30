from netCDF4 import Dataset
import numpy as np
import datetime as datetime
import calendar
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import stats
import xlwt,xlrd
from xlutils.copy import copy
import os

# ************************************************
#  Excel
# ************************************************

def output(ireg,reg,filename, sheet, title,month_abbr, cc):
    pad = 12+4

    if os.path.exists(filename):
        rb = xlrd.open_workbook(filename, formatting_info=True)
        sheets = rb.sheet_names()
        wb = copy(rb)
        if sheet in sheets:
            # r_sheet = rb.sheet_by_name(sheet)
            # print r_sheet
            sh = wb.get_sheet(sheets.index(sheet))
        else:
            sh = wb.add_sheet(sheet)
            for im in np.arange(1,12+1):
                sh.write(im+1,0,month_abbr[im])
                sh.write(pad+im,0,month_abbr[im])
            sh.write(0,0,title)
            sh.write(1,0,"corr")
            sh.write(pad,0,"corr_dt")

    else:
        wb = xlwt.Workbook()
        sh = wb.add_sheet(sheet)
        for im in np.arange(1,12+1):
            sh.write(im+1,0,month_abbr[im])
            sh.write(pad+im,0,month_abbr[im])
        sh.write(0,0,title)
        sh.write(1,0,"corr")
        sh.write(pad,0,"corr_dt")


    sh.write(1,ireg+1,reg)
    sh.write(pad,ireg+1,reg)

    for im in np.arange(1,12+1):
        # sh.write(im,ireg+0,month_abbr[im])
        sh.write(im+1,ireg+1,cc[0,im-1])
        sh.write(pad+im,ireg+1,cc[1,im-1])

    wb.save(filename)

def output_trend(ireg,reg,filename, sheet, title,month_abbr, trend):
    pad = 12+4

    if os.path.exists(filename):
        rb = xlrd.open_workbook(filename, formatting_info=True)
        sheets = rb.sheet_names()
        wb = copy(rb)
        if sheet in sheets:
            # r_sheet = rb.sheet_by_name(sheet)
            # print r_sheet
            sh = wb.get_sheet(sheets.index(sheet))
        else:
            sh = wb.add_sheet(sheet)
            for im in np.arange(1,12+1):
                sh.write(im+1,0,month_abbr[im])
                sh.write(pad+im,0,month_abbr[im])
            sh.write(0,0,title)
            sh.write(1,0,"slope")
            sh.write(pad,0,"p value")

    else:
        wb = xlwt.Workbook()
        sh = wb.add_sheet(sheet)
        for im in np.arange(1,12+1):
            sh.write(im+1,0,month_abbr[im])
            sh.write(pad+im,0,month_abbr[im])
        sh.write(0,0,title)
        sh.write(1,0,"slope")
        sh.write(pad,0,"p value")


    sh.write(1,ireg+1,reg)
    sh.write(pad,ireg+1,reg)

    for im in np.arange(1,12+1):
        # sh.write(im,ireg+0,month_abbr[im])
        sh.write(im+1,ireg+1,trend[0,im-1])
        # if trend[1,im-1] < .05:
        sh.write(pad+im,ireg+1,trend[1,im-1])

    wb.save(filename)



# ************************************************
#  Selection
# ************************************************
month_abbr = ["","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"]

# time
year = [1979,2016]
yrs = np.arange(year[0],year[1],1)
nyrs = np.size(yrs)
x = np.arange(0,nyrs-1,1)

mm = [1,12]  #[ 9,9]
dd = [1,31]
hr = [0,18]

# ssn   / use "YYY" for whole year
# ssn="YYY"

# nf
nf1 = 0# 5518 ;0
nf2 = 20000# 5518 ;20000

hs = "SH"
# reg
nreg = 5

nfig = 10  # possible combinations of vars
fig1 = 8
fig2 = 10
cc = np.zeros((nreg+1,nfig+1,2,12),dtype=np.float) # correlation


for ireg in np.arange(nreg+1):
    if ireg == 0:
        lon = [ -90,361]
        lat = [ -40,-20]
        reg = "20_40S"

    if ireg == 1:
        lon = [ 30,90]
        lat = [ -40,-20]
        reg = "30_90.20_40S"

    if ireg == 2:
        lon = [ 90,150]
        lat = [ -40,-20]
        reg = "90_150.20_40S"

    if ireg == 3:
        lon = [ 150,210]
        lat = [ -40,-20]
        reg = "150_210.20_40S"

    if ireg == 4:
        lon = [ 210,285]
        lat = [ -40,-20]
        reg = "210_285.20_40S"

    if ireg == 5:
        lon = [ 300,340]
        lat = [ -40,-20]
        reg = "300_340.20_40S"



    maxnf = 200 # max # of fronts per timestep
    maxnp = 100 # max # of frontal points


    # ************************************************
    # read in netCDF files
    # ************************************************
    # 1. slp
    # ************************************************

    for yr in yrs:
      slpyr = yr
      fin = "/Users/irudeva/work/DATA/ERAint/Mslp_highres/erain.mslp.avm.%d.nc"%slpyr
      print fin
      nc = Dataset(fin, 'r')

      # varnam=['longitude','latitude','time','msl']
      # v=0
      # for var in varnam:
      #   if nc.variables[varnam[v]].name != var:
      #       print "Variables don't agree", var, nc.variables[varnam[v]].name, v
      #       exit()
      #   v += 1

      # nc_attrs, nc_dims, nc_vars = ncdump(nc)

      # print nc_attrs, nc_dims, nc_vars
      print nc.variables.keys()
      if nc.variables.keys()[0]=="lon":
          lonslp = nc.variables['lon'][:]
          latslp = nc.variables['lat'][:]
      elif nc.variables.keys()[0]=="longitude":
          lonslp = nc.variables['longitude'][:]
          latslp = nc.variables['lattitude'][:]
      else:
         print"Check lon/lat names in the netcdf file, %d "%yr
      time = nc.variables['time'][:]
      mslp = nc.variables['msl'][:]/100.

      dt_nc = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
           for t in time]
      print dt_nc[1].year


      latslpSH = latslp[np.logical_and(latslp<-20,latslp>-65)]
      mslpSH   = mslp[:,np.logical_and(latslp<-20,latslp>-65),:]

      # zonal average
      slp_z = np.mean(mslpSH[:,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=2)

      if yr==year[0]:
          STR1slp = np.zeros((2,12,nyrs),dtype=np.float)
          STR1lat = np.zeros_like(STR1slp)
          STR1slp_tr = np.zeros((2,12),dtype=np.float)
          STR1lat_tr = np.zeros_like(STR1slp_tr)

          STR2slp = np.zeros_like(STR1slp)
          STR2lat = np.zeros_like(STR2slp)
          STR2slp_tr = np.zeros_like(STR1slp_tr)
          STR2lat_tr = np.zeros_like(STR1slp_tr)


      for im in np.arange(slp_z[:,0].size):
        #   print dt_nc[im].month
          STR1slp[0,im,yr-year[0]] = np.amax(slp_z[im,:])
          STR1lat[0,im,yr-year[0]] = latslpSH[np.argmax(slp_z[im,:])]
        #   print STR1lat[im,yr-year[0]], STR1slp[im,yr-year[0]]

          STR2slp[0,im,yr-year[0]] = np.mean(np.amax(mslpSH[im,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1))
          STR2lat[0,im,yr-year[0]] = np.mean(latslpSH[np.argmax(mslpSH[im,:,np.logical_and(lonslp<lon[1],lonslp>lon[0])],axis=1)])

      nc.close()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  Plotting
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # plt.close('all')
    #
    # f, ax = plt.subplots(6, 2,  sharex='col', sharey='row')
    # plt.suptitle('STR1 intensity, SH',fontsize=14)
    # for ic in np.arange(2):
    #     for ir in np.arange(6):
    #         im = ir+ic*6
    #         a = ax[ir, ic]
    #         a2 = a.twinx()
    #
    #         a.plot(yrs,STR1slp[im,:],color='b',lw=1,label='STR1 intensity')
    #         a2.plot(yrs,STR1lat[im,:],color='g',lw=1,label='STR1 location')
    #         # ax[ir, ic].set_title('STR1 intensity, %s, SH'%month_abbr[ir+ic*6+1])
    #         a.axis((year[0]-1, year[1], 1010, 1025))
    #         a2.axis((year[0]-1, year[1], -40, -25))
    #
    #         if im == 3:
    #             ax[ir, ic].set_ylabel('MSL Pressure (Pa)',color='b')
    #         if im == 9:
    #             a2.set_ylabel('latitude',color='g')
    #
    #         #  trend
    #         slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR1slp[im,:])
    #         trend = intercept + (slope * yrs)
    #         a.plot(yrs,trend,color='b',lw=2)
    #
    #         slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR1lat[im,:])
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

    # for yr in yrs:
    #   fin = "../frontgrd/frline.%d.nc"%yr
    #   print "fronts from %s"%fin
    #   ncf = Dataset(fin, 'r')
    #
    #   print ncf.variables.keys()
    #   flon = ncf.variables['flon'][:]
    #   flat = ncf.variables['flat'][:]
    #   time = ncf.variables['time'][:]
    #   dv = ncf.variables['dv'][:]
    #   nf = ncf.variables['nfr'][:]
    #   npts = ncf.variables['npts'][:]
    #
    #   fdt = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
    #        for t in time]
    #
    #
    #   if any(x > maxnf) :
    #     print "ERROR: nf > maxnf"
    #     quit()
    #
    #   if any(x > maxnp ) :
    #     print "ERROR:  np > maxnp"
    #     quit()
    #
    #   # ************************************************
    #   # masking
    #   # ************************************************
    #   # create mask
    #
    #   for ind,t in enumerate(fdt):
    #     if t.year == yr and t.month==1 and t.day == 1 and t.hour == 0:
    #        tfr1 = ind
    #     if t.year == yr and t.month==12 and t.day == 31 and t.hour == 18:
    #        tfr2 = ind
    #        break
    #
    #   if yr == year[0]:
    #       if calendar.isleap(yr):
    #           frmask = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
    #           frdv  = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
    #       else:
    #           frmask = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.int)
    #           frdv   = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.float32)
    #       frNPlat = np.zeros_like(frdv)
    #
    #
    #   for nt in np.arange(tfr1,tfr2+1):
    #       if fdt[nt].year>=year[0] and fdt[nt].year<=year[1]:
    #         #   for im in range( 1,13):
    #         #       if fdt[nt].month>=mm[0] and fdt[nt].month<=mm[1]:
    #         #         #   print im
    #                   for ifr in range(nf[nt]):
    #                       if ifr+1>=nf1 and ifr+1<=nf2:
    #                           for ip in range(npts[nt,ifr]):
    #                               if flat[nt,ifr,ip]>=lat[0] and flat[nt,ifr,ip]<=lat[1]:
    #                                   if flon[nt,ifr,ip]>=lon[0] and flon[nt,ifr,ip]<=lon[1]:
    #                                       frmask[yr-year[0],nt-tfr1,ifr] = fdt[nt].month
    #                                       frdv[yr-year[0],nt-tfr1,ifr] = np.mean(dv[nt,ifr,:npts[nt,ifr]])
    #                                       frNPlat[yr-year[0],nt-tfr1,ifr] = np.amax(flat[nt,ifr,:npts[nt,ifr]])
    #                                     #   print yr-year[0],nt-tfr1,ifr, "   ",im,month_abbr[im], frmask[yr-year[0],nt-tfr1,ifr]
    #                                     #   print np.sum(frmask==im)
    #                                       break
    #
    #
    #
    # # ************************************************
    # # front stats
    # # ************************************************
    #
    # dvcrit = np.zeros(12,dtype = np.float)
    # nfr_my    = np.zeros([ 2,12,nyrs],dtype=np.float)
    # nfr_my3tr = np.zeros_like(nfr_my)
    # frNPlat_my = np.zeros_like(nfr_my)
    # frNPlat_my3tr = np.zeros_like(nfr_my)
    #
    # nfr_my_tr   = np.zeros([ 2,12],dtype=np.float)
    # nfr_my3tr_tr = np.zeros_like(nfr_my_tr)
    # frNPlat_my_tr= np.zeros_like(nfr_my_tr)
    # frNPlat_my3tr_tr = np.zeros_like(nfr_my_tr)
    #
    # for im in range( 1,13):
    #
    #   frdvtmp = frdv[np.where(frmask == im)]
    #
    #   dvcrit[im-1] =  np.percentile(frdvtmp,67)
    #   print  month_abbr[im], dvcrit[im-1]
    #
    #   # hist, bin_edges = np.histogram(frdvtmp,range=(0.1,25))
    #   # print "pdf: ",hist
    #   # print bin_edges
    #
    #
    # for yr in yrs:
    #     for im in range( 1,13):
    #         nfr_my[0,im-1,yr-year[0]]=np.sum(frmask[yr-year[0],:,:]==im)
    #
    #         tmp0 = frNPlat[yr-year[0],:,:][np.where(frmask[yr-year[0],:,:] == im)]
    #         frNPlat_my[0,im-1,yr-year[0]]=np.mean(tmp0)
    #
    #         tmp  = frdv[yr-year[0],:,:][np.where(frmask[yr-year[0],:,:] == im)]
    #         nfr_my3tr[0,im-1,yr-year[0]]=np.count_nonzero(np.where(tmp>=dvcrit[im-1], 1,0))
    #         frNPlat_my3tr[0,im-1,yr-year[0]]=np.mean(tmp0[np.where(tmp>=dvcrit[im-1])])
    #         # nfr_my3tr[im-1,yr-year[0]]=np.sum(frmask[yr-year[0],:,:]==im and frdv[yr-year[0],:,:]>=dvcrit[im-1])
    #         # print yr, im, month_abbr[im], frNPlat_my[0,im-1,yr-year[0]],frNPlat_my3tr[0,im-1,yr-year[0]]
    #         # print np.count_nonzero(np.where(frmask[yr-year[0],:,:]== im, 1,0))

    # ************************************************
    # trends
    # ************************************************

    for im in range( 1,13):
        # slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,nfr_my[0,im-1,:])
        # nfr_my[1,im-1,:] = intercept + (slope * yrs)
        # nfr_my_tr[0,im-1] = slope
        # nfr_my_tr[1,im-1] = p_value
        # title = "Trends in N of fronts"
        # output_trend(ireg,reg,"../output/trend.Nfr.%d_%d.xls"%(year[0],year[1]-1), "trend",title, month_abbr, nfr_my_tr)
        #
        # slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,nfr_my3tr[0,im-1,:])
        # nfr_my3tr[1,im-1,:] = intercept + (slope * yrs)
        # nfr_m3try_tr[0,im-1] = slope
        # nfr_m3try_tr[1,im-1] = p_value
        # title = "Trends in N of STRONG fronts "
        # output_trend(ireg,reg,"../output/ntrend.Nfrstr.%d_%d.xls"%(year[0],year[1]-1), "trend", title,month_abbr, nfr_my3tr_tr)
        #
        # slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,frNPlat_my[0,im-1,:])
        # frNPlat_my[1,im-1,:] = intercept + (slope * yrs)
        # frNPlat_my_tr[0,im-1] = slope
        # frNPlat_my_tr[1,im-1] = p_value
        # title = "Trends in frontal NP lat"
        # output_trend(ireg,reg,"../output/trend.frNPlat.%d_%d.xls"%(year[0],year[1]-1), "trend", title,month_abbr, frNPlat_my_tr)
        #
        # slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,frNPlat_my3tr[0,im-1,:])
        # frNPlat_my3tr[1,im-1,:] = intercept + (slope * yrs)
        # frNPlat_my3tr_tr[0,im-1] = slope
        # frNPlat_my3tr_tr[1,im-1] = p_value
        # title = "Trends in STRONG fronts NP lat"
        # output_trend(ireg,reg,"../output/trend.frNPlat_str.%d_%d.xls"%(year[0],year[1]-1), "trend", title, month_abbr, frNPlat_my3tr_tr)

        slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR1lat[0,im-1,:])
        STR1lat[1,im-1,:] = intercept + (slope * yrs)
        STR1lat_tr[0,im-1] = slope
        STR1lat_tr[1,im-1] = p_value
        title = "Trends in STR1 location"
        output_trend(ireg,reg,"../output/trend.STR1loc.%d_%d.xls"%(year[0],year[1]-1), "trend", title,month_abbr, STR1lat_tr)

        slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR1slp[0,im-1,:])
        STR1slp[1,im-1,:] = intercept + (slope * yrs)
        STR1slp_tr[0,im-1] = slope
        STR1slp_tr[1,im-1] = p_value
        title = "Trends in STR1 intensity"
        output_trend(ireg,reg,"../output/trend.STR1int.%d_%d.xls"%(year[0],year[1]-1), "trend", title,month_abbr, STR1slp_tr)

        slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR2lat[0,im-1,:])
        STR2lat[1,im-1,:] = intercept + (slope * yrs)
        STR2lat_tr[0,im-1] = slope
        STR2lat_tr[1,im-1] = p_value
        title = "Trends in STR2 location"
        output_trend(ireg,reg,"../output/trend.STR2loc.%d_%d.xls"%(year[0],year[1]-1), "trend",title, month_abbr, STR2lat_tr)

        slope, intercept, r_value, p_value, std_err = stats.linregress(yrs,STR2slp[0,im-1,:])
        STR2slp[1,im-1,:] = intercept + (slope * yrs)
        STR2slp_tr[0,im-1] = slope
        STR2slp_tr[1,im-1] = p_value
        title = "Trends in STR2 intensity"
        output_trend(ireg,reg,"../output/trend.STR2int.%d_%d.xls"%(year[0],year[1]-1), "trend", title,month_abbr, STR2slp_tr)


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Plotting
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    majorLocator = MultipleLocator(5)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1)


    # plt.ion()
    plt.close('all')
    for ifig in np.arange(fig1,fig2+1):
        # fig = plt.figure(ifig)
        f, ax = plt.subplots(6, 2,  sharex='col', sharey='row',figsize=(3.5*3.13,4.2*3.13))

        if ifig ==0:
            title = 'STR intensity vs number of fronts, %s'%reg
            fout = "frN_STRint"
        if ifig ==1:
            title = 'STR intensity vs number of strong fronts, %s'%reg
            fout = "frNstr_STRint"
        if ifig ==2:
            title = 'STR location vs number of fronts, %s'%reg
            fout = "frN_STRloc"
        if ifig ==3:
            title = 'STR location vs number of strong fronts, %s'%reg
            fout = "frNstr_STRloc"
        if ifig ==4:
            title = 'STR intensity vs north lat of fronts, %s'%reg
            fout = "frNorthLat_STRint"
        if ifig ==5:
            title = 'STR intensity vs north lat of strong fronts, %s'%reg
            fout = "frNorthLatstr_STRint"
        if ifig ==6:
            title = 'STR location vs north lat of fronts, %s'%reg
            fout = "frNorthLat_STRloc"
        if ifig ==7:
            title = 'STR location vs north lat of strong fronts, %s'%reg
            fout = "frNorthLatstr_STRloc"
        if ifig ==8:
            title = 'STR1 intensity vs location, %s'%reg
            fout = "STR1int_STR1loc"
        if ifig ==9:
            title = 'STR1 intensity vs STR2 intensity, %s'%reg
            fout = "STR1int_STR2int"
        if ifig ==10:
            title = 'STR1 location vs STR2 location, %s'%reg
            fout = "STR1loc_STR2loc"
        plt.suptitle(title,fontsize=14)

        if ifig ==0 or ifig == 2:
            var1 = nfr_my
            yax1 = [900, 1800]
            yax1label = 'number of fronts'
            var1legend = 'N of fronts'

        if ifig ==1 or ifig == 3:
            var1 = nfr_my3tr
            yax1 = [300, 700]
            yax1label = 'number of fronts'
            var1legend = 'N of strong fronts (int > 67th perc)'

        if any ([4,6] == ifig ):
            var1 = frNPlat_my
            yax1 = [-40, -25]
            yax1label = 'latitude'
            var1legend = 'front lat north'

        if any ([5,7] == ifig ):
            var1 = frNPlat_my3tr
            yax1 = [-40, -25]
            yax1label = 'latitude'
            var1legend = 'strong front lat north'

        if any ([0,1,4,5] == ifig ):
            var2 = STR1slp
            yax2 = [1010, 1025]
            yax2label = 'MSL Pressure (Pa)'
            var2legend = 'STR intendity'

        if any ([0,1,4,5] == ifig ):
            var2 = STR1slp
            yax2 = [1010, 1025]
            yax2label = 'MSL Pressure (Pa)'
            var2legend = 'STR intendity'

        if any ([2,3,6,7,8] == ifig ):
            var2 = STR1lat
            yax2 = [-40, -25]
            yax2label = 'latitude'
            var2legend = 'STR location'

        if any ([8,9] == ifig ):
            var1 = STR1slp
            yax1 = [1010, 1025]
            yax1label = 'MSL Pressure (Pa)'
            var1legend = 'STR1 intensity'

        if ifig == 9:
            var2 = STR2slp
            yax2 = [1010, 1025]
            yax2label = 'MSL Pressure (Pa)'
            var2legend = 'STR2 intensity'

        if ifig == 10:
            var1 = STR1lat
            yax1 = [-40, -25]
            yax1label = 'latitude'
            var1legend = 'STR1 location'

            var2 = STR2lat
            yax2 = [-40, -25]
            yax2label = 'latitude'
            var2legend = 'STR2 location'


        for ic in np.arange(2):
            for ir in np.arange(6):
                im = ir+ic*6
                a = ax[ir, ic]
                a2 = a.twinx()

                a.plot(yrs,var1[0,im,:],color='b',lw=1,label=var1legend)
                a2.plot(yrs,var2[0,im,:],color='g',lw=1,label=var2legend)
                # ax[ir, ic].set_title('STR intensity, %s, SH'%month_abbr[ir+ic*6+1])

                # a.axis((year[0]-1, year[1], yax1[0], yax1[1]))
                # a2.axis((year[0]-1, year[1], yax2[0], yax2[1]))
                # ax2.set_ylim(0, 35)
                # ax.set_ylim(-20,100)

                # a.xaxis.set_ticks(yrs)
                # xa = a.get_xaxis()
                # a.xaxis.set_tick_params(axis='x',which='major',length=6,width =2)
                # a.xaxis.set_tick_params(axis='x',which='minor',length=3,width =2)

                a.xaxis.set_major_locator(majorLocator)
                a.xaxis.set_major_formatter(majorFormatter)
                a.yaxis.set_major_formatter(majorFormatter)

                # for the minor ticks, use no labels; default NullFormatter
                a.xaxis.set_minor_locator(minorLocator)

                a.grid()

                if im == 3:
                    ax[ir, ic].set_ylabel(yax1label,color='b')
                if im == 9:
                    a2.set_ylabel(yax2label,color='g')

                #  trend
                a.plot(yrs,var1[1,im,:],color='b',lw=2)
                a2.plot(yrs,var2[1,im,:],color='g',lw=2)

                cctmp = np.corrcoef(var1[0,im,:],var2[0,im,:])
                cc[ireg,ifig,0,im] = cctmp[1,0]
                cctmp = np.corrcoef(var1[0,im,:]-var1[1,im,:],var2[0,im,:]-var2[1,im,:])
                cc[ireg,ifig,1,im] = cctmp[1,0]

                a.set_title(' %s, %s, r = %.2f, r_dt = %.2f'%(month_abbr[ir+ic*6+1],reg,cc[ireg,ifig,0,im],cc[ireg,ifig,1,im]))

                if im==11:
                    # ask matplotlib for the plotted objects and their labels
                    lines, labels = a.get_legend_handles_labels()
                    lines2, labels2 = a2.get_legend_handles_labels()
                    a2.legend(lines + lines2, labels + labels2, loc='lower right')

                    # ax[ir, ic].legend(loc='lower right')

        plt.setp([a.set_xlabel('Year') for a in ax[5, :]])
        # a=ax[-1, -1]
        # ax[0, 0].xlabel('Year')
        # ax[0, 0].ylabel('MSL Pressure (Pa)')
        # ax[0, 0].axis((1979, 2014, 1012, 1030))


        f.subplots_adjust(hspace=0.3)
        plt.draw()
        f.savefig("../output/%s.%s.%d_%d.png"%(fout,reg,year[0],year[1]-1))
        plt.close(f)

        # ************************************************
        #  correlations to Excel
        # ************************************************
        fxls = "../output/%s.%d_%d.xls"%(fout,year[0],year[1]-1)
        if ireg == 0 and os.path.exists(fxls):
                os.remove(fxls)
        output(ireg,reg,fxls, "corr", title, month_abbr, cc[ireg,ifig,:,:])
