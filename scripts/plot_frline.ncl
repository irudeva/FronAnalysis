begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time
year = (/2015,2015/)
yr = 2015 ; year(0)
mm = (/9,9/)  ;(/9,9/)
dd = (/10,10/)
hr = (/18,18/)

;ssn   / use "YYY" for whole year
ssn="SON"

;nf
nf1 = 0; 5518 ;0
nf2 = 20000; 5518 ;20000

hs = "SH"
;reg
lon = (/10,150/)
lat = (/-60,-30/)
reg = "smth"


;************************************************
; read in netCDF files
;************************************************
; 1. slp
;************************************************

;t = 0
  ;slpyr = yr
  ;a = addfile("/Users/irudeva/work/DATA/ERAint/erain.mslp."+slpyr+".nc","r")

  ;;slp = short2flt(a->msl(t,:,:))
  ;slp = short2flt(a->msl(:,:,:))
  ;slp1 =0.01*slp
  ;copy_VarMeta(slp, slp1)

  ;timeslp = a->time
  ;dateslp = cd_calendar(timeslp, 0)

  ;print(dateslp)
  ;status_exit(1)

  ;do nt=0, dimsizes(timeslp)-1
  ;  if (dateslp(nt,0).eq.yr.and.dateslp(nt,1).eq.mm.and.dateslp(nt,2).eq.dd.and.dateslp(nt,3).eq.hr)then
  ;   tslp = nt
  ;   break
  ;  end if
  ;end do

  ;print(n)
  ;print (cd_calendar(timeslp(n), 0))
  ;status_exit(1)

  ;************************************************
  ; 2. fronts
  ;************************************************

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

;************************************************
; Plotting - preparation
;************************************************

;create mask
  frmask = flon
  frmask = 0

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
res@mpMinLatF            = -90      ; range to zoom in on
res@mpMaxLatF            =  -10.
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

;************************************************
; masking
;************************************************


do nt=0, dimsizes(time)-1
  ;print ("time="+datefr(nt,:))
  ;detect season
  if(any((/1,2,12/).eq.datefr(nt,1)))then
    tssn = "DJF"
  end if
  if(any((/3,4,5/).eq.datefr(nt,1)))then
    tssn = "MAM"
  end if
  if(any((/6,7,8/).eq.datefr(nt,1)))then
    tssn = "JJA"
  end if
  if(any((/9,10,11/).eq.datefr(nt,1)))then
    tssn = "SON"
  end if


  if (datefr(nt,0).ge.year(0).and.datefr(nt,0).le.year(1))then
   if (datefr(nt,1).ge.mm(0).and.datefr(nt,1).le.mm(1))then
    if (datefr(nt,2).ge.dd(0).and.datefr(nt,2).le.dd(1))then
     if (datefr(nt,3).ge.hr(0).and.datefr(nt,3).le.hr(1))then

      if (ssn.eq."YYY".or.tssn.eq.ssn)then

       do ifr = 0,nf(nt)-1
       ;print(ifr+" "+nf(nt))

        if (ifr+1.ge.nf1.and.ifr+1.le.nf2)then

         do np = 0,npts(nt,ifr)-1
          if (flat(nt,ifr,np).ge.lat(0).and.flat(nt,ifr,np).le.lat(1))then
           if (flon(nt,ifr,np).ge.lon(0).and.flon(nt,ifr,np).le.lon(1)) then
            frmask(nt,ifr,np) = 1
            ;print (nf+"  "+np)
            ;print (flon(nt,nf,n))
            ;print (flat(nt,nf,n))
            ;status_exit(1)
            end if  ;lon
           end if  ;lat

         end do   ; nf

        end if   ;fr1/fr2

       end do  ;

      end if  ; ssn

     end if ;hr
    end if ;dd
   end if ;mm
  end if ;year

end do ; time

;status_exit(1)


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
