begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time

 ;ssn   / use "YYY" for whole year
  ssn=(/"YYY","DJF","MAM","JJA","SON"/)

 year = (/1979,1983/)
 yr = year(0)
 mm = (/1,12/)  ;(/9,9/)
 dd = (/0,31/)
 hr = (/0,18/)

 if(mm(0).ne.1.and.ssn.eq."YYY")then
  print("WARNING: Check the season and time settings!")
  print("  season = "+ssn)
  print("  month(s) = "+month_abbr(mm(0))+" to "+month_abbr(mm(1)))
  print("  Check 'cyr' var!!!")
  sleep(1)
 end if

;ntrk
trk1 = 0; 5518 ;0
trk2 = 20000; 5518 ;20000

;reg

hs ="SH"

rlon = (/-181,361/)
rlat = (/-90,90/)
reg = ""


;rlon = (/100,150/)
;rlat = (/-60,-30/)
;reg = "AU"+"."

;************************************************
; Read STR loc
;************************************************

fSTR = "../output/ts/STRloc_5deg.20_40S.1979_2015.txt"
print ("read STRloc from "+fSTR)
z1 = readAsciiTable(fSTR,6,"float",1)

STRlon = z1(:,0)
; print(dimsizes(STRlon))
; print(dimsizes(STRlon(0)))
wSTRlon = array_append_record (STRlon, STRlon(0), 0)
; print (wSTRlon)
STRlat = new((/dimsizes(STRlon),dimsizes(ssn)/),float)
wSTRlat = new((/dimsizes(wSTRlon),dimsizes(ssn)/),float)

do issn=0,4
  ; print (z1(:,issn+1))
  STRlat(:,issn) = z1(:,issn+1)
  ; print (  STRlat(:,issn))
  wSTRlat(:,issn) = array_append_record (STRlat(:,issn), STRlat(0,issn), 0)

end do

; status_exit(1)
;************************************************
; Plotting - preparation
;************************************************

date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
             sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
             sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

regstr = reg+hs
wks = gsn_open_wks("png","../output/frnum."+"."+date_str1+"_"+date_str2+"."+regstr )                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame


res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = True       ; Turn on/off contour lines
res@cnFillPalette       = "WhiteBlueGreenYellowRed"   ; Set color map

res@gsnPolar   = hs

;res@mpLimitMode="LatLon"
if (hs.eq."NH")then
 res@mpMinLatF            = 25      ; range to zoom in on
 res@mpMaxLatF            =  90.
 else if (hs.eq."SH")then
   res@mpMinLatF            = -90      ; range to zoom in on
   res@mpMaxLatF            =  -20.
 else
   print("ERROR: Check hs")
 end if
 end if
;res@mpMinLonF           = -60.
;res@mpMaxLonF           = 120.


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF   = 1
res@cnMaxLevelValF   = 7
res@cnLevelSpacingF  = .5

res@lbLabelBarOn      = False           ; turn off individual cb's
res@cnLineLabelsOn    = False
res@cnLinesOn         = False                   ; turn off contour lines


; res@tiMainString        = date_str1 +" to "+ date_str2
; res@tiMainString        = "From "+date_str1+" to "+date_str2 + "  " + ssn+"~C~~Z90~ reg lon: ("+rlon(0)+"E to "+rlon(1)+"E)  lat: ("+rlat(0)+"N to "+rlat(1)+"N)"

plres = True
plres@gsLineColor = "red"
plres@gsLineThicknessF = 3.0

plot = new(5,graphic)
plotSTR = new(5,graphic)

;************************************************
; read in netCDF files
;************************************************

do issn = 1,3,2

  do iyr =year(0),year(1)

  frin = "../frontgrd/frgrd."+iyr+".nc"
  print("reading file = "+frin)
  b = addfile(frin,"r")
  frloc = b->fronts(:,:,:)
  frarea = b->wfronts(:,:,:)   ;the same as frlocbut 3 grid point wide, for smoother plots
  lat = b->latitude
  lon = b->longitude


  timecyc = b->time
  datecyc = cd_calendar(timecyc, 0)

  ;print(datecyc)
  ;status_exit(1)

  if(iyr.eq.year(0))then

   frgrd = new((/dimsizes(lat),dimsizes(lon)/),float)
   frgrd = 0
   frgrd!1 = "lon"
   frgrd!0 = "lat"

   frgrd&lon = lon
   frgrd&lat = lat

   ntstp = 0

  end if

  do nt=0, dimsizes(timecyc)-1
    cyr = datecyc(nt,0)
    if(any((/1,2,12/).eq.datecyc(nt,1)))then
      if(datecyc(nt,1).eq.12.and.ssn(issn).ne."YYY") then
       cyr = cyr+1
      end if
      cssn = "DJF"
    end if
    if(any((/3,4,5/).eq.datecyc(nt,1)))then
      cssn = "MAM"
    end if
    if(any((/6,7,8/).eq.datecyc(nt,1)))then
      cssn = "JJA"
    end if
    if(any((/9,10,11/).eq.datecyc(nt,1)))then
      cssn = "SON"
    end if
  ;print(datecyc(nt,:))

    if (cyr.ge.year(0).and.cyr.le.year(1).and.cyr.eq.iyr)then
     if (datecyc(nt,1).ge.mm(0).and.datecyc(nt,1).le.mm(1))then
      if (datecyc(nt,2).ge.dd(0).and.datecyc(nt,2).le.dd(1))then
       if (datecyc(nt,3).ge.hr(0).and.datecyc(nt,3).le.hr(1))then

        if (ssn(issn).eq."YYY".or.cssn.eq.ssn(issn))then

          frgrd = frgrd+where(frarea(nt,:,:).gt.0, 1, 0)
          ntstp  = ntstp + 1

         end if

        end if ;hr
       end if ;dd
      end if ;mm
     end if ;year
    ;end do  ;nit
  ;  end if ;
  ;  print ("ku-ku")
  end do

  delete(frloc)
  delete(frarea)
  delete(timecyc)
  delete(datecyc)

  end do  ; year

  ;************************************************
  ; 5. Plotting
  ;************************************************

 print ("Start plotting")

 res@tiMainString        = ssn(issn)

 ; twist!!!
  ; tmpgrd = frgrd
  ; tmpgrd(:,:299) = frgrd(:,60:)
  ; tmpgrd(:,300:) = frgrd(:,:59)
  ; frgrd = tmpgrd
  ; frgrd&lon = lon
 ;end twist


 ;tstr = cd_calendar(timecyc(tcyc), 3)

 ;month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
 ;                     "Oct","Nov","Dec"/)
 ;year   = tointeger(datecyc(tcyc,0))    ; Convert to integer for
 ;month  = tointeger(datecyc(tcyc,1))    ; use sprinti
 ;day    = tointeger(datecyc(tcyc,2))
 ;hour   = tointeger(datecyc(tcyc,3))
 ;date_strm = sprinti("%0.4i", year)+month_abbr(month)+ \
 ;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)
 ;date_str = sprinti("%0.4i", year)+sprinti("%0.2i", month)+ \
;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)

 ;print(date_str + "")
 ;status_exit(1)

frgrd = frgrd*100/ntstp

regfrgrd = frgrd
regfrgrd = 0
regfrgrd(ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1))) = frgrd(ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1)))

wrf_smooth_2d(regfrgrd,3)
plot(issn) = gsn_csm_contour_map_polar(wks,regfrgrd,res)
;g = smth9(regfrgrd, 0.50, -0.25, True)
;g = smth9(g, 0.50, -0.25, True)
;copy_VarMeta(frgrd,g)
; plot1 = gsn_csm_contour_map_polar(wks,g,res)


plotSTR(issn) = gsn_add_polyline(wks, plot(issn), wSTRlon, wSTRlat(:,issn), plres)
;  draw(plot(issn))
; frame(wks)

 end do  ; issn
;************************************************
; create panel
;************************************************
  resP                    = True                 ; modify the panel plot
  ; resP@gsnPanelMainString = "A common title"     ; new resource added in NCL V6.4.0
  ;resP@txString           = "fr N: "+STR +" "+sel(0) +"-"+ sel(1)+", "+ ssn + "  "+ regtitle
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  resP@lbLabelFontHeightF  = 0.007               ; make labels smaller
  resP@lbLabelStride       = 2


  gsn_panel(wks,plot(1:3:2),(/1,2/),resP)               ; now draw as one plot

delete(plot)



status_exit(1)

end
