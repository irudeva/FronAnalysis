from netCDF4 import Dataset
import numpy as np
# from windspharm.standard import VectorWind
# from windspharm.tools import prep_data, recover_data, order_latdim
# from scipy import interpolate
import datetime as datetime
from datetime import timedelta  # Python standard library datetime  module
# import time as ftime
# import scipy.ndimage.filters as filters
# import scipy.ndimage.morphology as morphology
# import scipy.ndimage as sp
import calendar


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def dist(plat,plon,lat,lon):

    rad  =6371032 #in m
    pi   = np.pi
    dtr  = pi/180.
    rtd  = 180./pi

    dist=0.
    lat1=plat*dtr
    lon1=plon*dtr

    lat2=lat*dtr
    lon2=lon*dtr

    if lon1==lon2 :
        dist=abs(lat1-lat2)
        dist=rad*dist
    elif lat1==lat2 :
        dist=np.arccos(np.cos(lat1)**2*(np.cos(lon2-lon1)-1.)+1.)
        dist=rad*dist
    else :
        dist=np.arccos(np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
        dist=rad*dist

    return dist


##  time - t-1 ???
## what is in the last column of fronts?

R    = 6371032 #in m
pi   = np.pi
deg  = R*pi/180

#create lon/lat grid
grdres = 1  # grid reolution
lons = np.arange(0,360,grdres)
lats = np.arange(90,-91,-grdres)

frwidth = 2 # number of point that are used to 'inflate' fronts to the west

dset = "erain"
# scheme = "UOM"
hs   = ["0N","0S"]

fyr = 2015
lyr = 2015

yr0 = 1979 # the first year in the dataset

#IMPORTANT: for trk files from 16 Nov YR-1 to 31 Jan YR+1
# netcdf reference time
time0 = 'hours since 1900-01-01 00:00:0.0'
dt0 = datetime.datetime(1900, 1, 1, 0)

max_nf      =  250     # max number of tracks per year
max_npts    =  200      # max number of points per front
min_flength =  500      # min front length
inf_npts    =  2        #number of points to inflate fronts to the west

# creat arrays for recording tracks
#trklen  = np.empty(max_ntrk,np.int16); trklen.fill(0)
#trktime = np.empty([max_ntrk,max_trklength],np.int); trktime.fill(0)
#trklon  = np.empty([max_ntrk,max_trklength]); trklon.fill(0)
#trklat  =  np.copy(trklon); trklat.fill(0)
#trkslp  =  np.copy(trklon); trklat.fill(0)
#trkrd   =  np.copy(trklon); trklat.fill(0)
#trklpl  =  np.copy(trklon); trklat.fill(0)
#trkdp   =  np.copy(trklon); trklat.fill(0)
#trkiop  = np.copy(trktime); trklat.fill(0)


for yr in range(fyr,lyr+1):
    #create time variable
    if yr == yr0 :
        dt_nc = [ datetime.datetime(yr, 1, 1, 0)+i*timedelta(hours=6) for i in range(1460)]
        if calendar.isleap(yr):
            dt_nc = [ datetime.datetime(yr, 1, 1, 0)+i*timedelta(hours=6) for i in range(1464)]
    else :
        dt_nc = [ datetime.datetime(yr-1, 11, 1, 0)+i*timedelta(hours=6) for i in range(1764)]
        if calendar.isleap(yr):
            dt_nc = [ datetime.datetime(yr-1, 11, 1, 0)+i*timedelta(hours=6) for i in range(1768)]

    dt_hours = np.empty(len(dt_nc)); dt_hours.fill(0)
    for it, dt in enumerate(dt_nc) :
        tdelta =  dt - datetime.datetime(1900, 1, 1, 0)
        dt_hours[it] = divmod(tdelta.total_seconds(),3600)[0]


    #grid for front records
    nf_dt  = np.empty(len(dt_nc),np.int); nf_dt.fill(0)  # number of fronts per time step
    frgrd = np.empty([len(lons),len(lats),len(dt_nc)],np.int16); frgrd.fill(0)
    inffrgrd = np.copy(frgrd)
    # cycgrd_rad  = np.copy(frgrd)    #for radius from cyc tracking
    dvgrd  = np.copy(frgrd,np.float)

    cglon  = np.zeros([len(dt_nc),max_nf])
    cglat  = np.zeros_like(cglon)
    clen  = np.zeros_like(cglon)
    npts  = np.zeros_like(cglon,dtype = np.int16)

    klon = np.zeros([len(dt_nc),max_nf,max_npts])
    klat = np.zeros_like(klon)
    kdv  = np.zeros_like(klon)


    #read front
    for ihs in hs :
        fin = "../../FrontTRK/output/erain/10/fcyc2trk10.%s.%d.%s.dat"%(dset,yr,ihs)
        print "fronts from =",fin
        print "Reading fronts:"

        f = open(fin, 'r')

        for tind,t in enumerate(dt_nc) :
        # for line in range(1,max_ntrk):
            emptyline = f.readline()
            header    = f.readline()
            emptyline = f.readline()
            smthline = f.readline()
            emptyline = f.readline()

            # fout.write(emptyline)
            # fout.write(header)
            # fout.write(emptyline)
            # fout.write(smthline)

            if header == '':
             print ' fin is at the eof'
             break
            else :
             print ' header: ', header.strip()

            columns = header.split()

            #check the date
            cdate=columns[1]
            cyr = int(cdate[:2])
            cm  = int(cdate[2:4])
            cd  = int(cdate[4:6])
            chr = int(columns[2][:2])
            if cyr < 30:
                if any(cyr == iyr for iyr in range(18,23)) & any(yr == iyr for iyr in range(1998,2002)) :
                    cyr = cyr + 1980
                else :
                    cyr = cyr + 2000
            else:
                cyr = cyr + 1900

            if t != datetime.datetime(cyr,cm,cd,chr):
                print "ERROR: Dates do not agree"
                print " t = ",t
                print columns[1:3]
                print datetime.datetime(cyr,cm,cd,chr)
                quit()

            # N of fronts
            cnf = int(columns[6])

            # done to here!
            nf_dt[tind] = nf_dt[tind]+cnf   #combined number of fronts (for two henispheres)
            nf = nf_dt[tind]
            if nf > max_nf :
                print " ERROR!!! Number of fronts is more than max_nf: nf = %d, max_nf = %d" %(nf,max_nf)
                print "          time = ", t
                quit()

            for ifr in range(nf-cnf,nf):
                # print "t: %s ifr=%d / %d" % (t,ifr+1,nf)
                l = f.readline()
                columns = l.split()
                if int(columns[0])-1  != ifr-(nf-cnf) :
                    print "ERROR: Check the front number"
                    print "       ifr = ",int(columns[0]), ", current number = ", ifr-(nf-cnf)+1
                    print "       time: ", t
                    quit()
                cglon[tind,ifr] = float(columns[2])
                cglat[tind,ifr] = float(columns[3])
                clen[tind,ifr]  = float(columns[4])
                npts[tind,ifr]  = float(columns[7])

                if npts[tind,ifr] > max_npts:
                    print "ERROR: the front %d is too long" % ifr
                    print "       t: ",t
                    quit()

                for k in range(npts[tind,ifr]) :
                    iline    = f.readline()
                    columns = iline.split()

                    klon[tind,ifr,k] = float(columns[4])
                    klat[tind,ifr,k] = float(columns[3])
                    kdv[tind,ifr,k]  = float(columns[5])


                ####################################################
                # fornts onto grid

                #find the closest grid value to the cyclone center
                # clon[ntrk-1,n] = 1 # test
                # clat[ntrk-1,n] = 55.  #test!
                if clen[tind,ifr] >= min_flength :

                    # interpolation front lines onto grid-associated lon/lat
                    cpts = npts[tind,ifr]
                    y1 = find_nearest(lats,klat[tind,ifr,0])
                    y2 = find_nearest(lats,klat[tind,ifr,cpts-1])

                    newcpts = abs(y2-y1)+1

                    # check 0 line crossing
                    if any(abs(klon[tind,ifr,i+1]-klon[tind,ifr,i])>300 for i in range(cpts-1)):
                        print  "ERROR: front potentially crosses the zero line, check interpolation !"
                        print "       time : ", t, " ifr = ", ifr
                        print "       lons : ",klon[tind,ifr,:cpts]
                        print "       diff : ",[abs(klon[tind,ifr,i+1]-klon[tind,ifr,i]) for i in range(cpts-1)]
                        quit()

                    if lats[y1]>lats[y2] :
                        ynew = np.arange(lats[y1],lats[y2]-grdres,-grdres)
                        xnew = np.interp(ynew, klat[tind,ifr,cpts-1::-1], klon[tind,ifr,cpts-1::-1])
                    elif lats[y1]<lats[y2] :
                        ynew = np.arange(lats[y1],lats[y2]+grdres,grdres)
                        xnew = np.interp(ynew, klat[tind,ifr,:cpts], klon[tind,ifr,:cpts])
                    else:
                        print "ERROR: check y1 and y2"
                        print "       time : ", t, " ifr = ", ifr
                        print "        lat : ",klat[tind,ifr,:cpts]
                        print "        y1 = %d, y2 = %d" % (y1,y2)
                        quit()


                    # for k in range(npts[tind,ifr]) :
                    #
                    #     if klon[tind,ifr,k] < 0:
                    #         klon[tind,ifr,k] = klon[tind,ifr,k]+360
                    #     ilon = find_nearest(lons,klon[tind,ifr,k])
                    #     ilat = find_nearest(lats,klat[tind,ifr,k])
                    #
                    #     frgrd[ilon,ilat,tind] = ifr-(nf-cnf)+1   #ifr+1
                    #     dvgrd[ilon,ilat,tind]  = kdv[tind,ifr,k]

                    for k in range(newcpts) :

                        if xnew[k] < 0:
                            # print "--", tind,ifr,k,xnew[k]
                            xnew[k] = xnew[k]+360
                            # print "!+",xnew[k],ynew[k]
                            # print find_nearest(lons, xnew[k] )
                            # print lons[find_nearest(lons, xnew[k])]
                        ilon = find_nearest(lons, xnew[k] )
                        ilat = find_nearest(lats, ynew[k] )

                        frgrd[ilon,ilat,tind] = ifr-(nf-cnf)+1   #ifr+1
                        dvgrd[ilon,ilat,tind]  = kdv[tind,ifr,k]


                        # inflation
                        if inf_npts == 2 :
                            if ilon > 1 :
                                inffrgrd[ilon-2:ilon+1,ilat,tind] = ifr+1
                            elif ilon > 0 :
                                inffrgrd[ilon-1:ilon+1,ilat,tind] = ifr+1
                                inffrgrd[-1,ilat,tind] = ifr+1
                            elif ilon == 0 :
                                inffrgrd[-2:,ilat,tind] = ifr+1
                                # print "inf: ",inffrgrd[-2:,ilat,tind]
                        else:
                            print "ERROR: Check the number of points to inflate fronts"


            # break   #for the first time step only
        f.close()
        print 'fin closed'

    #---Cyclone number to NetCDF write---------------------------------------------------------------
    print("Start NetCDF writing")



    # ncvar = "cycnum"
    # print 'ncvar=',ncvar
    frfile = '../frontgrd/frgrd.%d.nc' % (yr)
    print "    write to %s " % frfile
    ncfr = Dataset(frfile, 'w', format='NETCDF4')
    ncfr.description = "Front lines from  %s. Min front length is %d " % (fin,min_flength)

    dimnam=('lon','lat','time','nf')
    varnam=['longitude','latitude','time','nf','fronts','wfronts',"dv"]

    #dimensions
    ncfr.createDimension(dimnam[0], lons.size)
    ncfr.createDimension(dimnam[1], lats.size)
    ncfr.createDimension(dimnam[2], None)
    ncfr.createDimension(dimnam[3], len(dt_nc))


    #variables
    # for nv in range(0, 3) :
    #     # ncfr_var = ncfr.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
    #     ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[nv])
        # for ncattr in nc.variables[varnam[nv]].ncattrs():
            # ncfr_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
    #print(nc.variables['latitude'].ncattrs())

    #create  variables
     #lons
    nv=0
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[nv])
    if varnam[nv] == 'longitude' :
        ncfr_var.long_name = varnam[nv]
        ncfr_var.units = 'degrees_east'
        ncfr.variables[varnam[nv]][:] = lons

     #lats
    nv=1
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[nv])
    if varnam[nv] == 'latitude' :
        ncfr_var.long_name = varnam[nv]
        ncfr_var.units = 'degrees_north'
        ncfr.variables[varnam[nv]][:] = lats

     #time
    nv=2
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[nv])
    ncfr_var.long_name = varnam[nv]
    if varnam[nv] == 'time' :
        ncfr_var.calendar = 'gregorian'
        # ncfr_var.units = 'hours since 1900-01-01 00:00:0.0'
        ncfr_var.units = time0
        #for one time step - test!
        # tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
        # ncfr.variables[varnam[nv]][:] = divmod(tdelta.total_seconds(),3600)[0]
        ncfr.variables[varnam[nv]][:] = dt_hours

    #front to netcdf
    nv = 3
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[2])
    # print ncfr_var.shape
    ncfr_var.long_name = 'number of fronts'
    # ncfr_var.scale_factor = varlist["scale"][iv]
    # ncfr_var.add_offset   = 0.
    # ncfr_var.units        = 'scale   %s' % varlist["scale"][iv]

    ncfr_var[:] = nf_dt


    nv = 4
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[-2::-1])
    # print ncfr_var.shape
    ncfr_var.long_name = 'location of fronts'
    # ncfr_var.scale_factor = varlist["scale"][iv]
    # ncfr_var.add_offset   = 0.
    # ncfr_var.units        = 'scale   %s' % varlist["scale"][iv]

    #print qx.shape
    #print ncfr_var.shape
    # ncfr_var[:,:,:] = np.swapaxes(gridcyc,0,2)
    ncfr_var[:,:,:] = np.swapaxes(frgrd[:,:,:],0,2)

    #cinflated fronts to netcdf
    nv =5
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[-2::-1])
    # print ncfr_var.shape
    ncfr_var.long_name = 'inflated fronts'
    ncfr_var.width = "plus %d grid points to the west" % inf_npts
    # ncfr_var.scale_factor = varlist["scale"][iv]
    # ncfr_var.add_offset   = 0.
    # ncfr_var.units        = 'scale   %s' % varlist["scale"][iv]

    ncfr_var[:,:,:] = np.swapaxes(inffrgrd[:,:,:],0,2)

    nv = 6
    ncfr_var = ncfr.createVariable(varnam[nv], 'f',dimnam[-2::-1])
    # print ncfr_var.shape
    ncfr_var.long_name = 'wind shift'
    # ncfr_var.scale_factor = varlist["scale"][iv]
    # ncfr_var.add_offset   = 0.
    ncfr_var.units        = 'm/s'

    ncfr_var[:,:,:] = np.swapaxes(frgrd[:,:,:],0,2)




    ncfr.close()

    ##---Tracks to NetCDF write---------------------------------------------------------------


    frlines = '../frontgrd/frline.%d.nc' % (yr)
    print "    write to %s " % frlines
    ncfline = Dataset(frlines, 'w', format='NETCDF4')
    ncfline.description = "Front lines from  %s" % (fin)

    #dimensions
    dimnam=('n','nf','time')

    ncfline.createDimension(dimnam[0], max_npts )
    ncfline.createDimension(dimnam[1], max_nf )
    ncfline.createDimension(dimnam[2], None)

    #variables
    #  1D variable
    ncvar = ncfline.createVariable('time', 'f',dimnam[2])
    ncvar.calendar = 'gregorian'
    # ncfr_var.units = 'hours since 1900-01-01 00:00:0.0'
    ncvar.units = time0
    ncvar[:] = dt_hours


    ncvar = ncfline.createVariable('nfr', 'i4',dimnam[2])
    ncvar.long_name = 'number of fronts per timestep'
    ncvar[:]    = nf_dt

    #  2D variables
    ncvar = ncfline.createVariable('npts', 'i2',dimnam[:0:-1])
    ncvar.long_name = 'number of points per front'
    ncvar[:,:]    = npts



    frmt = "(%s,%s,%s)f4" % (len(dt_nc),max_nf,max_npts)
    varlist = np.zeros( 3, dtype = {  'names': ['name', 'long_name', 'dtype',  'units', 'data'],
                                    'formats': [  'a7',       'a31',    'a5',    'a44',   frmt]})


    varlist[0] = ( "flon","longitude of frontal points",  'f', 'degrees_east', klon[:,:,:])
    varlist[1] = ( "flat", "latitude of frontal points",  'f','degrees_north', klat[:,:,:])
    varlist[2] = (   "dv",                 "wind shift",  'f',          'm/s',  kdv[:,:,:])


    # varnam=['trklen','trktime','trklon','trklat','trkslp','trkrad','trklpl','trkdp','trkiop']

    for iv in range(varlist['name'].size) :
        # print varlist["name"][iv]
        ncvar = ncfline.createVariable(varlist["name"][iv], varlist["dtype"][iv],dimnam[::-1])
        ncvar.long_name = varlist["long_name"][iv]
        ncvar.units     = varlist[    "units"][iv]
        ncvar[:]        = varlist[     "data"][iv]


    ncfline.close()


    ##---End NetCDF write---------------------------------------------------------------
