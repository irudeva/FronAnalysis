begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time
year = (/1979,2015/)
yrs = ispan(year(0),year(1),1)
nyrs = dimsizes(yrs)
x = ispan(0,nyrs-1,1)

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
reg = "20_40S"

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

  timefr= b->time
  datefr = cd_calendar(timefr, 0)

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

  do nt=0, dimsizes(timefr)-1
    if (datefr(nt,0).eq.yr.and.datefr(nt,1).eq.1.and.datefr(nt,2).eq.1.and.datefr(nt,3).eq.0)then
       tfr1 = nt
    end if
    if (datefr(nt,0).eq.yr.and.datefr(nt,1).eq.12.and.datefr(nt,2).eq.31.and.datefr(nt,3).eq.18)then
       tfr2 = nt
       break
    end if
  end do

  if(yr .eq. year(0))then
   if(isleapyear(yr))then
    frmask = new((/dimsizes(yrs),tfr2-ntfr1+1,maxnf,maxnp/),integer)
    frdv  = new((/dimsizes(yrs),tfr2-tfr1+1,maxnf/),float)
   else
    frmask = new((/dimsizes(yrs),tfr2-tfr1+5,maxnf,maxnp/),integer)
    ; frdv  = new((/dimsizes(yrs),tfr2-tfr1+5,maxnf/),float)
    frdv = new((/dimsizes(yrs),tfr2-tfr1+5,maxnf,maxnp/),float)
   end if
   frmask = 0
   frdv  = 0
 end if


 ; do nt=0, dimsizes(timefr)-1
 do nt=tfr1,tfr2

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
              ;print(yr+" "+year(0)+" "+nt+" "+ifr+" "+np)

              ; for all fronts that had at least one point in the region
              frmask(yr-year(0),nt-tfr1,ifr,np) = im
              frdv (yr-year(0),nt-tfr1,ifr,np) = avg(dv(nt,ifr,:npts(nt,ifr)-1))

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

  delete(timefr)
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

  frmask1d =  ndtooned(frmask)
  frdv1d  =  ndtooned(frdv)

  nfr_my    = new((/2,12,nyrs/),float)  ; the first dim is added for trends
  nfr_my3tr = nfr_my

  do im =1,12

   print(month_abbr(im)+"")
   ; print (im+" "+ind(frmask1d.eq.im))
   nn = num(frmask1d.eq.im)
   zpdf = pdfx(frdv1d(ind(frmask1d.eq.im)),0,False)
   ; printVarSummary( zpdf )

   print(zpdf@bin_center+" "+zpdf*100/nn)

   opt = False
   stat = stat_dispersion(frdv1d(ind(frmask1d.eq.im)), opt )
   ; for high tritile
   print ("high tritile = "+stat(9))

  do yr = year(0),year(1)

    frmaskyr = ndtooned(frmask(yr-year(0),:,:,:))
    frdvyr = ndtooned(frdv(yr-year(0),:,:,:))

    nfr_my(0,im-1,yr-year(0)) = num(frmaskyr.eq.im)
    ; for high tritile
    frmasktmp = where(frdvyr.ge.stat(9),frmaskyr,0)
    nfr_my3tr(0,im-1,yr-year(0)) = num(frmasktmp.eq.im)
    ;print(yr+" "+nfr_my(im-1,yr-year(0)) + "  "+num(frmaskyr.eq.im))

    ; calculated trends
    rc    = regline(x,nfr_my(0,im-1,:))          ; slope
    rc@units = "front/month"
    nfr_my(1,im-1,:) = rc*x + rc@yintercept

    rc    = regline(x,nfr_my3tr(0,im-1,:))          ; slope
    rc@units = "front/month"
    nfr_my3tr(1,im-1,:) = rc*x + rc@yintercept


  end do  ;yr
  end do  ;im
;************************************************
; Plotting
;************************************************

plot = new(12,graphic)

date_str1 = sprinti("%0.4i", year(0))   ;+sprinti("%0.2i", mm(0)))    ;+ \
          ;   sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))   ;+sprinti("%0.2i", mm(1)))    ;+ \
          ;   sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

wks = gsn_open_wks("png","../output/frN_STR."+date_str1+"_"+date_str2+"."+reg)                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

res@tmXBMode          = "Explicit"              ; explicit labels
res@tmXBValues        =  x               ; location of labels
res@tmXBLabels        =  yrs              ; labels themselves
res@tmLabelAutoStride = True                    ; nice stride on labels

res@xyMarkLineModes     = (/"Lines","Lines"/)  ; choose which have markers
;res@xyMarkers           = 16                     ; choose type of marker
;res@xyMarkerColor       = "red"                  ; Marker color
res@xyLineColor         = "red"                  ; Marker color
;res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)
;res@xyDashPatterns      = 1                      ; solid line
res@xyLineThicknesses   = (/1,2/)                ; set second line to 2
;res@tmYLFormat          = "f"                    ; not necessary but nicer labels




do im =0,11
  res@tiMainString      = month_abbr(im+1)
  plot(im)  = gsn_csm_xy (wks,x,nfr_my3tr(:,im,:),res) ; create plot
end do

;************************************************
; create panel
;************************************************
  resP                    = True                 ; modify the panel plot
  ; resP@gsnPanelMainString = "A common title"     ; new resource added in NCL V6.4.0
 resP@txString           = "Number of strong fronts vs STR lat from "+date_str1+" to "+date_str2 + "  "
  gsn_panel(wks,plot,(/6,2/),resP)               ; now draw as one plot


end
