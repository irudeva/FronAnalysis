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
year = [1979,1984]#2016]
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
lon = [ -90,361]
lat = [ -40,-20]
reg = "20_40S"

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
#       STRslp = np.zeros((12,nyrs),dtype=np.float)
#       STRlat = np.zeros_like(STRslp)
#
#   for im in np.arange(slp_z[:,0].size):
#     #   print dt_nc[im].month
#       STRslp[im,yr-year[0]] = np.amax(slp_z[im,:])
#       STRlat[im,yr-year[0]] = latslpSH[np.argmax(slp_z[im,:])]
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
  fin = "../frontgrd/frline.%d.nc"%yr
  print "fronts from %s"%fin
  ncf = Dataset(fin, 'r')

  print ncf.variables.keys()
  flon = ncf.variables['flon'][:]
  flat = ncf.variables['flat'][:]
  time = ncf.variables['time'][:]
  dv = ncf.variables['dv'][:]
  nf = ncf.variables['nfr'][:]
  npts = ncf.variables['npts'][:]

  fdt = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
       for t in time]


  if any(x > maxnf) :
    print "ERROR: nf > maxnf"
    quit()

  if any(x > maxnp ) :
    print "ERROR:  np > maxnp"
    quit()

  # ************************************************
  # masking
  # ************************************************
  # create mask

  for ind,t in enumerate(fdt):
    if t.year == yr and t.month==1 and t.day == 1 and t.hour == 0:
       tfr1 = ind
    if t.year == yr and t.month==12 and t.day == 31 and t.hour == 18:
       tfr2 = ind
       break

  if yr == year[0]:
      if calendar.isleap(yr):
          frmask = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.int)
          frdv  = np.zeros((nyrs,tfr2-tfr1+1,maxnf),dtype=np.float)
      else:
          frmask = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.int)
          frdv   = np.zeros((nyrs,tfr2-tfr1+5,maxnf),dtype=np.float32)


  for nt in np.arange(tfr1,tfr2+1):
      if fdt[nt].year>=year[0] and fdt[nt].year<=year[1]:
        #   for im in range( 1,13):
        #       if fdt[nt].month>=mm[0] and fdt[nt].month<=mm[1]:
        #         #   print im
                  for ifr in range(nf[nt]):
                      if ifr+1>=nf1 and ifr+1<=nf2:
                          for ip in range(npts[nt,ifr]):
                              if flat[nt,ifr,ip]>=lat[0] and flat[nt,ifr,ip]<=lat[1]:
                                  if flon[nt,ifr,ip]>=lon[0] and flon[nt,ifr,ip]<=lon[1]:
                                      frmask[yr-year[0],nt-tfr1,ifr] = fdt[nt].month
                                      frdv[yr-year[0],nt-tfr1,ifr] = np.mean(dv[nt,ifr,:npts[nt,ifr]])
                                    #   print yr-year[0],nt-tfr1,ifr, "   ",im,month_abbr[im], frmask[yr-year[0],nt-tfr1,ifr]
                                    #   print np.sum(frmask==im)
                                      break



# ************************************************
# front stats
# ************************************************

dvcrit = np.zeros(12,dtype = np.float)
nfr_my    = np.zeros([ 12,nyrs],dtype=np.float)
nfr_my3tr = np.zeros_like(nfr_my)


for im in range( 1,13):

  frdvtmp = frdv[np.where(frmask == im)]

  dvcrit[im-1] =  np.percentile(frdvtmp,67)
  print  month_abbr[im], dvcrit[im-1]

  # hist, bin_edges = np.histogram(frdvtmp,range=(0.1,25))
  # print "pdf: ",hist
  # print bin_edges


for yr in yrs:
    for im in range( 1,13):
        nfr_my[im-1,yr-year[0]]=np.sum(frmask[yr-year[0],:,:]==im)
        tmp  = frdv[yr-year[0],:,:][np.where(frmask[yr-year[0],:,:] == im)]
        nfr_my3tr[im-1,yr-year[0]]=np.count_nonzero(np.where(tmp>=dvcrit[im-1], 1,0))
        # nfr_my3tr[im-1,yr-year[0]]=np.sum(frmask[yr-year[0],:,:]==im and frdv[yr-year[0],:,:]>=dvcrit[im-1])
        print yr, im, month_abbr[im], nfr_my[im-1,yr-year[0]],nfr_my3tr[im-1,yr-year[0]]
        # print np.count_nonzero(np.where(frmask[yr-year[0],:,:]== im, 1,0))

# # ************************************************
# # Plotting
# # ************************************************
#
# plot = new(12,graphic)
# data = new([ 4,nyrs],float)
#
# date_str1 = sprinti("%0.4i", year(0))   ;+sprinti("%0.2i", mm(0)))    ;+ \
#           #   sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
# date_str2 = sprinti("%0.4i", year(1))   ;+sprinti("%0.2i", mm(1)))    ;+ \
#           #   sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))
#
# wks = gsn_open_wks("png","../output/frNstrong_STR2."+date_str1+"_"+date_str2+"."+reg)                  # send graphics to PNG file
#
# resF = True
# resF@gsnDraw              = False             # do not draw the plot
# resF@gsnFrame             = False             # do not advance the frame
#
# resF@tmXBMode          = "Explicit"              # explicit labels
# resF@tmXBValues        =  x               # location of labels
# resF@tmXBLabels        =  yrs              # labels themselves
# resF@tmLabelAutoStride = True                    # nice stride on labels
#
# resF@xyMarkLineModes     = [ "Lines","Lines"]  # choose which have markers
# ;resF@xyMarkers           = 16                     # choose type of marker
# ;resF@xyMarkerColor       = "red"                  # Marker color
# resF@xyLineColor         = "red"                  # Marker color
# ;resF@xyMarkerSizeF       = 0.005                  # Marker size (default 0.01)
# ;resF@xyDashPatterns      = 1                      # solid line
# resF@xyLineThicknesses   = [ 1,2]                # set second line to 2
# ;resF@tmYLFormat          = "f"                    # not necessary but nicer labels
#
# resP = True
# resP = resF
#
# resP@xyLineColor         = "blue"                  # Marker color
#
# res_text               = True
# res_text@txFontHeightF = 0.03                       # change font size
#
# amres                  = True
# amres@amJust           = "BottomCenter"
# amres@amParallelPosF   =  0.0    # This is the center of the plot.
# amres@amOrthogonalPosF = -0.72   # This is above the top edge of the plot.
#
#
#
#
# do im =0,11
#   resF@tiMainString      = month_abbr(im+1)
#   plot(im)  = gsn_csm_xy2 (wks,x,nfr_my3tr(:,im,:),STRlat(:,:,im),resF,resP) # create plot
#   # plot(im)  = gsn_csm_xy2 (wks,x,nfr_my(:,im,:),STRlat(:,:,im),resF,resP) # create plot
#
#   # correlation and significance
#   r    = escorc(nfr_my3tr(0,im,:),STRlat(0,:,im))
#   # r    = escorc(nfr_my(0,im,:),STRlat(0,:,im))
#   t    = r*sqrt((nyrs-2)/(1-r^2))
#   p    = student_t(t, nyrs-2)
#   psig = 0.05                       # test significance level
#   if (p.le.psig) then
#       text = "r="+r+" is significant at the 95% level"
#   else
#       text = "r="+r+" is NOT significant at the 95% level"
#   end if
#   # text_plot = gsn_create_text(wks, text, res_text)
#   print (month_abbr(im+1)+"  "+text)
#   ;gsn_add_annotation(plot(im), text_plot, amres)
#
#   # getvalues plot(im)
#   #  "tmYLLabelFontHeightF"   : fheight
#   #  "tmXTValues"             : tmXTValues
#   #  "tmYLValues"             : tmYLValues
#   # end getvalues
#   # nTm  = dimsizes(tmXTValues)               # number of major tick marks
#   # gsn_add_text(wks,bot_plot(im),text_plot,0.75*tmXTValues(nTm-1), \
#   #                                         0.35*tmYLValues(nTm-1) ,res_text)
#   gsn_add_text(wks,plot(im),text,1980,450 ,res_text)
#
#
# end do
#
# # ************************************************
# # create panel
# # ************************************************
#   print("Panel plot")
#   resall                    = True                 # modify the panel plot
#   # resP@gsnPanelMainString = "A common title"     # new resource added in NCL V6.4.0
#   resP@txString           = "Number of fronts vs STR lat from "+date_str1+" to "+date_str2 + "  "
#   gsn_panel(wks,plot,[ 6,2],resall)               # now draw as one plot
