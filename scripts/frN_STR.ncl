begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time
year = (/1979,1980/)
yrs = ispan(year(0),year(1),1)
mm = (/1,12/)  ;(/9,9/)
dd = (/1,31/)
hr = (/0,18/)

;ssn   / use "YYY" for whole year
;ssn="YYY"

;nf
nf1 = 0; 5518 ;0
nf2 = 20000; 5518 ;20000

hs = "SH"
;reg
lon = (/-90,361/)
lat = (/-40,-20/)
reg = "20_40"

maxnf = 200 ; max # of fronts per timestep
maxnp = 100 ; max # of frontal points


;************************************************
; read in netCDF files
;************************************************
; 1. slp
;************************************************

; do yr = year(0),year(1)
;   slpyr = yr
;   a = addfile("/Users/irudeva/work/DATA/ERAint/Mslp_1_5deg/erain.slp."+slpyr+".nc","r")
;
;   ;slp = short2flt(a->msl(t,:,:))
;   slp = short2flt(a->msl(:,:,:))
;   slp1 =0.01*slp
;   ; copy_VarMeta(slp, slp1)
;
;   timeslp = a->time
;   dateslp = cd_calendar(timeslp, 0)
;
;   vNames = getfilevarnames(a)
;
;   if (vNames(0).eq."lon")then
;      lonslp = a->lon
;      latslp = a->lat
;    else if (vNames(0).eq."longitude")then
;      lonslp = a->longitude
;      latslp = a->latitude
;    else
;      print("Check lon/lat names in the netcdf file, "+yr)
;    end if
;    end if
;
;   ; print(dateslp)
;   ; status_exit(1)
;
;   do ilat = 0,dimsizes(lat)
;     avslp = dim_avg_n_Wrap(slp1(:,:,ind(lonslp.ge.lon(0).and.lonslp.le.lon(1))),2)
;   end do
;
;   ; do im = 1,12
;   ;   do nt=0, dimsizes(timeslp)-1
;   ;    if (dateslp(nt,0).eq.yr.and.dateslp(nt,1).eq.mm)then
;   ;     avslp(im,:) =
;   ;     break
;   ;    end if
;   ;   end do
;   ; end do
;
;   if (yr .eq. year(0)) then
;     slp_zm = new((/12,dimsizes(latslp)/),float)
;     STRlat = new((/dimsizes(yrs),12/),float)
;     STRint = STRlat
;   end if
;
;   do im = 1,12
;     ; print(im+" ")
;     ;print(ind(dateslp(:,1).eq.im))
;     slp_zm(im-1,:) = dim_avg_n_Wrap(avslp(ind(dateslp(:,1).eq.im),:),0)
;     ; print(slp_zm(im-1,ind(latslp.le.0)))
;     ; print(max(slp_zm(im-1,ind(latslp.le.0))))
;     ; print(maxind(slp_zm(im-1,ind(latslp.le.0))))
;     STRlat(yr-year(0),im-1) = latslp(maxind(slp_zm(im-1,ind(latslp.le.0))))
;     STRint(yr-year(0),im-1) = max(slp_zm(im-1,ind(latslp.le.0)))
;     ; print (STRlat(yr-year(0),im-1)+"")
;     ; print (STRint(yr-year(0),im-1)+"")
;     print(yr +" " + month_abbr(im)+" "+STRlat(yr-year(0),im-1)+" "+STRint(yr-year(0),im-1))
;   end do
;
;   delete(slp)
;   delete(slp1)
;   delete(avslp)
;   delete(timeslp)
;   delete(dateslp)
;
; end do
; ; status_exit(1)

  ;************************************************
  ; 2. fronts
  ;************************************************

do yr = year(0),year(1)
  fin = "../frontgrd/frline."+yr+".nc"
  print ("fronts from "+fin)
  b = addfile(fin,"r")

  time= b->time
  datefr = cd_calendar(time, 0)

  flon = b->flon
  flat = b->flat
  dv = b->dv
  npts = b->npts
  nf = b->nfr
  ;status_exit(1)

  if(any(nf.gt.maxnf))then
    print("ERROR: nf > maxnf")
    status_exit(1)
  end if

  if(any(npts.gt.maxnp))then
    print("ERROR: np > maxnp")
    status_exit(1)
  end if



  ;************************************************
  ; masking
  ;************************************************
  ;create mask
  if(yr .eq. year(0))then
   frmask = new((/dimsizes(yrs),dimsizes(time)+4,maxnf,maxnp/),integer)
   frmask = 0
   frdv  = new((/dimsizes(yrs),dimsizes(time)+4,maxnf/),float)
   frdv  = 0
 end if


  do nt=0, dimsizes(time)-1

    ;print ("time="+datefr(nt,:))
    ;detect season

    ; if(any((/1,2,12/).eq.datefr(nt,1)))then
    ;   tssn = "DJF"
    ; end if
    ; if(any((/3,4,5/).eq.datefr(nt,1)))then
    ;   tssn = "MAM"
    ; end if
    ; if(any((/6,7,8/).eq.datefr(nt,1)))then
    ;   tssn = "JJA"
    ; end if
    ; if(any((/9,10,11/).eq.datefr(nt,1)))then
    ;   tssn = "SON"
    ; end if
    ;

    if (datefr(nt,0).ge.year(0).and.datefr(nt,0).le.year(1))then
    ;  if (datefr(nt,1).ge.mm(0).and.datefr(nt,1).le.mm(1))then
    ;   if (datefr(nt,2).ge.dd(0).and.datefr(nt,2).le.dd(1))then
    ;    if (datefr(nt,3).ge.hr(0).and.datefr(nt,3).le.hr(1))then

        ; if (ssn.eq."YYY".or.tssn.eq.ssn)then
      do im = 1,12
       if (datefr(nt,1).ge.im.and.datefr(nt,1).le.mm(1))then
         do ifr = 0,nf(nt)-1
           ;print(ifr+" "+nf(nt))

          if (ifr+1.ge.nf1.and.ifr+1.le.nf2)then

           do np = 0,npts(nt,ifr)-1
            if (flat(nt,ifr,np).ge.lat(0).and.flat(nt,ifr,np).le.lat(1))then
             if (flon(nt,ifr,np).ge.lon(0).and.flon(nt,ifr,np).le.lon(1)) then
              frmask(yr-year(0),nt,ifr,np) = im
              frdv (yr-year(0),nt,ifr) = avg(dv(nt,ifr,:npts(nt,ifr)-1))

              ;print (nf+"  "+np)
              ;print (flon(nt,nf,n))
              ;print (flat(nt,nf,n))
              ;status_exit(1)
              end if  ;lon
             end if  ;lat

           end do   ; nf

          end if   ;fr1/fr2

         end do  ; nf
       end if
      end do

        ; end if  ; ssn

    ;    end if ;hr
    ;   end if ;dd
    ;  end if ;mm
    end if ;year

  end do ; time
  delete(time)
  delete(datefr)
  delete(flon)
  delete(flat)
  delete(npts)
  delete(nf)
  delete(dv)

end do  ;years

  ;status_exit(1)
  ;************************************************
  ; front stats
  ;************************************************

  do im =1,12

  zpdf = pdfx(frdv(ind(frmask.eq.im)),0,False)
  printVarSummary( zpdf )

  ; do nt=0, dimsizes(time)-1
  ;   ; print (" plotting nt = "+nt+" "+dimsizes(time)+" "+datefr(nt,0)+month_abbr(tointeger(datefr(nt,1)))+datefr(nt,2)+" "+datefr(nt,3)+"h")
  ;   do ifr = 0,nf(nt)-1
  ;     if(any(frmask(yr-year(0),nt,ifr,:).eq.im))then
  ;       frnum(yr-year(0),im,)
  ;     end if
  ;   end do
  ; end do
  end do ; im
;************************************************
; Plotting - preparation
;************************************************



date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
             sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
             sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

wks = gsn_open_wks("png","../output/frline."+ssn+"."+date_str1+"_"+date_str2+"."+reg+"."+hs )                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

;res@cnFillOn            = True        ; Turn on contour fill
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

res@tiMainString        = "From "+date_str1+" to "+date_str2 + "~C~~Z90~ reg lon: ("+lon(0)+"E to "+lon(1)+"E)  lat: ("+lat(0)+"N to "+lat(1)+"N)"
res@tiMainFontHeightF   = 0.015

;plot1 = gsn_csm_contour_map(wks,slp1,res)
;plot1 = gsn_csm_contour_map_polar(wks,slp1(tslp,:,:),res)
plot1 = gsn_csm_map(wks,res)


;draw(plot1)
;frame(wks)


;status_exit(1)

plres = True

 plres@gsLineColor      = "red"                  ; color of lines
 plres@gsLineThicknessF = 2.0                    ; thickness of lines

mres = True

 mres@gsMarkerIndex = 16
 mres@gsMarkerSizeF = .003

txres               = True
 txres@txFontHeightF = 0.01

;
;************************************************
; Plotting
;************************************************

do nt=0, dimsizes(time)-1
  print (" plotting nt = "+nt+" "+dimsizes(time)+" "+datefr(nt,0)+month_abbr(tointeger(datefr(nt,1)))+datefr(nt,2)+" "+datefr(nt,3)+"h")
  do ifr = 0,nf(nt)-1
  if(any(frmask(nt,ifr,:).eq.1))then
    nit = npts(nt,ifr)

  ;************************************************
  ; add polyline
  ;************************************************

  ;print (ntrk)
  ;print (flon(ntrk,:nit-1))
  ;print (flat(ntrk,:nit-1))
  str = unique_string("poly")            ;-- result is poly0-poly15
  plot1@$str$ = gsn_add_polyline(wks,plot1,flon(nt,ifr,:nit-1),flat(nt,ifr,:nit-1),plres)

  ;************************************************
  ; add polymarker
  ;************************************************
  mres@gsMarkerColor      = 20
  str = unique_string("marker")            ;-- result is poly0-poly15
  ;plot1@$str$ = gsn_add_polymarker(wks,plot1,flon(nt,ifr,:nit-1),flat(nt,ifr,:nit-1),mres)

  ;************************************************
  ; add date
  ;************************************************

  it = 0  ;which date to show

  yy     = tointeger(datefr(nt,0))    ; Convert to integer for
  month  = tointeger(datefr(nt,1))    ; use sprinti
  day    = tointeger(datefr(nt,2))
  hour   = tointeger(datefr(nt,3))
  date2plot = sprinti("%0.4i", yy)+month_abbr(month)+ \
             sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)

  txt = unique_string("txt")
  ;plot1@$txt$ = gsn_add_text(wks,plot1,date2plot,flon(nt,ifr,it)+3.,flat(nt,ifr,it)-3.,txres)

  ;status_exit(1)
  ;break
 end if
 end do  ;ifront

end do  ;nt

draw(plot1)
frame(wks)

end
