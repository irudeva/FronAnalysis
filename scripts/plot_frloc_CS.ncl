begin
;************************************************
; read in netCDF files
;************************************************
; 1. slp
;************************************************
yr = 1979
mm = 9
dd = 7
hr = 18
;t = 0
hs ="NH"
  a = addfile("/Users/irudeva/work/DATA/ERAint/erain.mslp."+yr+".nc","r")

  ;slp = short2flt(a->msl(t,:,:))
  slp = short2flt(a->msl(:,:,:))
  slp1 =0.01*slp
  copy_VarMeta(slp, slp1)

  timeslp = a->time
  dateslp = cd_calendar(timeslp, 0)

  ;print(dateslp)
  ;status_exit(1)

  do nt=0, dimsizes(timeslp)-1
    if (dateslp(nt,0).eq.yr.and.dateslp(nt,1).eq.mm.and.dateslp(nt,2).eq.dd.and.dateslp(nt,3).eq.hr)then
     tslp = nt
     break
    end if
  end do

  ;print(n)
  ;print (cd_calendar(timeslp(n), 0))
  ;status_exit(1)

  ;************************************************
  ; 2. fronts
  ;************************************************
  if (mm.ge.10)then
     smon = ""+mm
  else
     smon = "0"+mm
  end if
  b = addfile("~/work/Communication/ClemensSpensberger/frontalzones/ei.ans."+yr+smon+".850.frozo_obj.nc","r")
  frloc = b->frozo_obj(:,0,:,:)
  lat = b->latitude
  lon = b->longitude

  ;printVarSummary(frloc)
  ;printVarSummary(slp1)
  ;printVarSummary(lat)
  ;copy_VarMeta(slp, frloc)

  timefr = b->time
  datefr = cd_calendar(timefr, 0)

  ;print(datefr)
  ;status_exit(1)

  do nt=0, dimsizes(timefr)-1
    if (datefr(nt,0).eq.yr.and.datefr(nt,1).eq.mm.and.datefr(nt,2).eq.dd.and.datefr(nt,3).eq.hr)then
     tfr = nt
     break
    end if
  end do
  ;print(n)
  ;print (cd_calendar(timefr(tfr), 0))
  ;status_exit(1)


 tstr = cd_calendar(timefr(tfr), 3)

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                     "Oct","Nov","Dec"/)
 year   = tointeger(datefr(tfr,0))    ; Convert to integer for
 month  = tointeger(datefr(tfr,1))    ; use sprinti
 day    = tointeger(datefr(tfr,2))
 hour   = tointeger(datefr(tfr,3))
 date_strm = sprinti("%0.4i", year)+month_abbr(month)+ \
            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)
 date_str = sprinti("%0.4i", year)+sprinti("%0.2i", month)+ \
            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)


print(date_str + "")
;status_exit(1)

;************************************************
; Cyclones
;************************************************

c = addfile("~/work/Projects/CycTrkUOM/ERAint/trkgrid/cycloc_rd."+yr+".nc","r")
cycloc = c->cycloc(:,:,:)
cycarea = c->cycarea(:,:,:)
clat = c->latitude
clon = c->longitude

;printVarSummary(cycloc)
;printVarSummary(slp1)
;printVarSummary(lat)
;copy_VarMeta(slp, cycloc)

timecyc = c->time
datecyc = cd_calendar(timecyc, 0)

;print(datecyc)
;status_exit(1)

do nt=0, dimsizes(timecyc)-1
  if (datecyc(nt,0).eq.yr.and.datecyc(nt,1).eq.mm.and.datecyc(nt,2).eq.dd.and.datecyc(nt,3).eq.hr)then
   tcyc = nt
   break
  end if
end do


;************************************************
; Plotting
;************************************************
wks = gsn_open_wks("png","../output/frloc_CS."+date_str+"."+hs )                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = True       ; Turn on/off contour lines
res@cnFillPalette       = "BlueDarkRed18"  ; Set color map

res@gsnPolar   = hs

;res@mpLimitMode="LatLon"
if (hs.eq."NH")then
 res@mpMinLatF            = 25      ; range to zoom in on
 res@mpMaxLatF            =  90.
 else if (hs.eq."SH")then
   res@mpMinLatF            = -90      ; range to zoom in on
   res@mpMaxLatF            =  -25.
 else
   print("ERROR: Check hs")
 end if
 end if

;res@mpMinLonF            =  -60.
;res@mpMaxLonF            = 120.


res@tiMainString        = date_strm

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 920
res@cnMaxLevelValF = 1050
res@cnLevelSpacingF = 10.

;plot1 = gsn_csm_contour_map(wks,slp1,res)
plot1 = gsn_csm_contour_map_polar(wks,slp1(tslp,:,:),res)
;plot1 = gsn_csm_map(wks,res)


;draw(plot1)
;frame(wks)
;************************************************
;Plot cyclones
;************************************************

mres = True

mres@gsMarkerIndex = 16
mres@gsMarkerSizeF = .003

txres               = True
txres@txFontHeightF = 0.01

do j = 0,dimsizes(cycloc(tcyc,:,0))-1
;print ("j="+j)

 indloc = ind(cycloc(tcyc,j,:).gt.0)
 indarea = ind(cycarea(tcyc,j,:).gt.0)
 ;print (indices)

 ;print(ismissing(indices))
 if(.not.any(ismissing(indarea)))then
  mlon = clon(indarea)
  mlat = new(dimsizes(mlon),float)
  mlat = clat(j)

  ;print (j+" mlon="+mlon)
  ;print (j+" mlat="+mlat)

  mres@gsMarkerColor      = "blue"

  str = unique_string("poly")            ;-- result is poly0-poly15
  plot1@$str$ = gsn_add_polymarker(wks,plot1,mlon,mlat,mres)

  if(.not.any(ismissing(indloc)))then

   cyclon = clon(indloc)
   cyclat = new(dimsizes(cyclon),float)
   cyclat = clat(j)
   ;print (j+" cyclon="+cyclon)
   ;print (j+" cyclat="+cyclat)


   mres@gsMarkerColor   = "red"
   str = unique_string("poly")            ;-- result is poly0-poly15
   plot1@$str$ = gsn_add_polymarker(wks,plot1,cyclon,cyclat,mres)
   ;print("cyclon"+cyclon)
   ;print("cyclat"+cyclat)
   ;gsn_polymarker(wks,plot1,cyclon,cyclat,mres)

   ; labling
   trknum = new(dimsizes(cyclon),"float")
  ;  print(dimsizes(cyclon))
  ;  print(cycloc(tcyc,j,indloc))
  ;  print(tcyc+" "+j+" "+indloc)
   trknum = cycloc(tcyc,j,indloc)

   ;print(trknum)


   txt = unique_string("txt")
   plot1@$txt$ = gsn_add_text(wks,plot1,trknum,cyclon+3.,cyclat-3.,txres)


   delete(cyclon)
   delete(cyclat)
   delete(trknum)

  end if

  delete(mlon)
  delete(mlat)
 end if
 delete(indarea)
 delete(indloc)

end do

;************************************************
; Plot fronts
;************************************************

mres = True

mres@gsMarkerIndex = 16
mres@gsMarkerSizeF = .003

;mres@gsnDraw              = False             ; do not draw the plot
;mres@gsnFrame             = False             ; do not advance the frame

;mres@cnFillOn            = False        ; Turn on contour fill
;mres@cnLinesOn           = True       ; Turn off contour lines
;mres@cnMonoFillColor     = True         ; Use single pattern color;
;mres@cnLineColor         = "green"      ; Turn off contour lines

txres               = True
txres@txFontHeightF = 0.01


do j = 0,dimsizes(frloc(tfr,:,0))-1
;print ("j="+j+"  num="+num(frloc(tfr,j,:).gt.0))

if(num(frloc(tfr,j,:).gt.0).gt.0) then

   ;indloc = ind(frloc(tfr,j,:).gt.0)
   indarea = ind(frloc(tfr,j,:).gt.0)

   ;to plot inflated cyclones
   if(.not.any(ismissing(indarea)))then
    mlon = lon(indarea)
    mlat = new(dimsizes(mlon),float)
    mlat = lat(j)

    ;print (j+" mlon="+mlon)
    ;print (j+" mlat="+mlat)

    mres@gsMarkerColor      = "green"

    str = unique_string("poly")            ;-- result is poly0-poly15
    plot1@$str$ = gsn_add_polymarker(wks,plot1,mlon,mlat,mres)


    ;to plot fronts as they are
    ; if(.not.any(ismissing(indloc)))then
    ;
    ;  wlon = lon(indloc)
    ;  wlat = new(dimsizes(wlon),float)
    ;  wlat = lat(j)
    ;
    ;  mres@gsMarkerColor   = 200
    ;  str = unique_string("poly")            ;-- result is poly0-poly15
    ;  ;print (j+" wlon="+wlon)
    ;  ;print (j+" wlat="+wlat)
    ;
    ;  plot1@$str$ = gsn_add_polymarker(wks,plot1,wlon,wlat,mres)
    ;
    ;  ; labling
    ;  frnum = new(dimsizes(wlon),"float")
    ;  ;print(dimsizes(mlon))
    ;  ;print(frloc(tfr,j,indloc))
    ;  ;print(tfr+" "+j+" "+indloc)
    ;  frnum = frloc(tfr,j,indloc)
    ;
    ;  ;print(trknum)
    ;
    ;
    ;  txt = unique_string("txt")
    ;  ;plot1@$txt$ = gsn_add_text(wks,plot1,frnum,mlon+3.,mlat-3.,txres)
    ;  delete(frnum)
    ;
    ;
    ;  delete(wlon)
    ;  delete(wlat)
    ;
    ; end if

    delete(mlon)
    delete(mlat)
   end if
   delete(indarea)
  ;  delete(indloc)
  end if   ; num frloc>0

end do

;contour lines for climatologies
;res1=True

;res1@gsnDraw              = False             ; do not draw the plot
;res1@gsnFrame             = False             ; do not advance the frame

;res1@cnFillOn            = False        ; Turn on contour fill
;res1@cnLinesOn           = True       ; Turn off contour lines
;res1@cnMonoFillColor     = True         ; Use single pattern color
;res1@cnLineColor         = "green"      ; Turn off contour lines

;res1@cnLevelSelectionMode = "ManualLevels"
;res1@cnMinLevelValF = 1.
;res1@cnMaxLevelValF = 23.
;res1@cnLevelSpacingF = 2.

;plot2 = gsn_csm_contour(wks,frloc,res1)         ; create a default plot


;overlay(plot1,plot2)


draw(plot1)
frame(wks)



end
