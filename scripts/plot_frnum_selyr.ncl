begin
;************************************************
; Estimates trends in frontal points
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time

;ssn   / use "YYY" for whole year
ssnarr=(/"YYY","DJF","MAM","JJA","SON"/)

hs =(/"NH","SH"/)
;hs = "NH"
do ihs =1,1

STRarr = (/"STRloc", "STRint"/)
sel = (/"max","min"/)
do istr = 0,1
  STR = STRarr(istr)
  do ireg=0,5
      if( ireg .eq. 0) then
          rlon = (/ -90,361/)
          rlat = (/ -40,-20/)
          reg = abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = hs(ihs)
      end if

      if (ireg .eq. 1)then
          rlon = (/ 30,90/)
          rlat = (/ -40,-20/)
          reg = abs(rlon(0))+"_"+abs(rlon(1))+"E."+abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = abs(rlon(0))+"_"+abs(rlon(1))+"E."
        end if

      if (ireg .eq. 2)then
          rlon = (/ 90,150/)
          rlat = (/ -40,-20/)
          reg = abs(rlon(0))+"_"+abs(rlon(1))+"E."+abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = abs(rlon(0))+"_"+abs(rlon(1))+"E."
        end if

      if (ireg .eq. 3)then
          rlon = (/ 150,210/)
          rlat = (/ -40,-20/)
          reg = abs(rlon(0))+"_"+abs(rlon(1))+"E."+abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = abs(rlon(0))+"_"+abs(rlon(1))+"E."
        end if

      if (ireg .eq. 4)then
          rlon = (/ 210,285/)
          rlat = (/ -40,-20/)
          reg = abs(rlon(0))+"_"+abs(rlon(1))+"E."+abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = abs(rlon(0))+"_"+abs(rlon(1))+"E."
        end if

      if (ireg .eq. 5)then
          rlon = (/ 300,340/)
          rlat = (/ -40,-20/)
          reg = abs(rlon(0))+"_"+abs(rlon(1))+"E."+abs(rlat(1))+"_"+abs(rlat(0))+"S"
          regtitle = abs(rlon(0))+"_"+abs(rlon(1))+"E."
        end if


  do issn = 0,4

    ssn = ssnarr(issn)


  do syr = 0,0
    if (syr.eq.0)then
     year = (/1979,2015/)
    end if
    if (syr.eq.1)then
      year = (/1992,2002/)
    end if
    if (syr.eq.2)then
      year = (/2005,2015/)
    end if

    ; read years
    fyrs = "../output/STRyrs/"+STR+"."+reg+"."+ssn+year(0)+"_"+year(1)+".txt"
    print ("read years from "+fyrs)
    z1 = asciiread(fyrs,(/10,4/),"float")

    yrs1 = z1(:,1)
    yrs2 = z1(:,3)


   ; yrs = ispan (year(0),year(1),1)
   ; yr = year(0)
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
  maxfrN = 200 ; max number of fronts per timestep




  ;************************************************
  ; Plotting - preparation
  ;************************************************

  ; date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
  ;              sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
  ; date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
  ;              sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

  regstr = reg+hs(ihs)
  ; wks = gsn_open_wks("png","../output/FNPnum_"+STR+"."+reg+"."+ssn+year(0)+"_"+year(1) )                  ; send graphics to PNG file
  wks = gsn_open_wks("png","../output/STRyrs/frNall_STR/frN_"+STR+"10."+reg+"."+ssn+year(0)+"_"+year(1) )                  ; send graphics to PNG file

  res = True
  resd = True
  res@gsnDraw              = False             ; do not draw the plot
  res@gsnFrame             = False             ; do not advance the frame


  res@cnFillOn            = True        ; Turn on contour fill
  res@cnFillPalette       = "WhiteBlueGreenYellowRed"  ; Set color map

  res@gsnPolar   = hs(ihs)

  ;res@mpLimitMode="LatLon"
  ; if (hs(ihs).eq."NH")then
  ;  res@mpMinLatF            = 15      ; range to zoom in on
  ;  res@mpMaxLatF            =  90.
  ;  else if (hs(ihs).eq."SH")then
  ;    res@mpMinLatF            = -90      ; range to zoom in on
  ;    res@mpMaxLatF            =  -15.
  ;  else
  ;    print("ERROR: Check hs")
  ;  end if
  ;  end if
  ;res@mpMinLonF           = -60.
  ;res@mpMaxLonF           = 120.

  res@cnLinesOn       = False                   ; turn off contour lines
  res@cnLineLabelsOn  = False                   ; tuen off line labels

  ; res@tiMainString        = "Number of fronts per "+ssn
  ;res@tiMainString        = date_str1 +" to "+ date_str2
  ;res@tiMainString        = "From "+date_str1+" to "+date_str2 + "  " + ssn+"~C~~Z90~ reg lon: ("+rlon(0)+"E to "+rlon(1)+"E)  lat: ("+rlat(0)+"N to "+rlat(1)+"N)"



  resd@cnLevelSelectionMode = "ManualLevels"
  if (ssn.eq."YYY")then
    res@cnMinLevelValF   = 0.
    res@cnMaxLevelValF   = 25.
    res@cnLevelSpacingF  = 2.
  else
    res@cnMinLevelValF   = 0.
    res@cnMaxLevelValF   = 8.
    res@cnLevelSpacingF  = .5
  end if

  ; for diff
  resd = res
  resd@cnFillPalette       = "BlueDarkRed18"  ; Set color map

  ;res@cnLevelSelectionMode = "AutomaticLevels"
  resd@cnLevelSelectionMode = "ManualLevels"
  if (ssn.eq."YYY")then
    resd@cnMinLevelValF   = -5.
    resd@cnMaxLevelValF   = 5.
    resd@cnLevelSpacingF  = .5
  else
    resd@cnMinLevelValF   = -2.
    resd@cnMaxLevelValF   = 2.
    resd@cnLevelSpacingF  = .2
  end if


  plot = new(3,graphic)                          ; create a plot array



  ;************************************************
  ; read in netCDF files
  ;************************************************

    do nsel =0,1
      if(nsel .eq. 0) then
         yrs = yrs1
      end if
      if(nsel .eq. 1) then
         yrs = yrs2
      end if
      ntstp = 0

      do iyr =0, dimsizes(yrs)-1

        yr = yrs(iyr)
        frin = "../frontgrd/frgrd."+yr+".nc"
        print("reading file = "+frin)
        b = addfile(frin,"r")
        frloc = b->fronts(:,:,:)
        frarea = b->wfronts(:,:,:)   ;the same as frlocbut 3 grid point wide, for smoother plots
        lat = b->latitude
        lon = b->longitude


        timent = b->time
        datent = cd_calendar(timent, 0)

        ;print(datent)
        ;status_exit(1)

        if(yr.eq.yrs(0).and.nsel.eq.0)then

         frgrd = new((/2,dimsizes(yrs),dimsizes(lat),dimsizes(lon)/),float)
         frgrd = 0
         frgrd!3 = "lon"
         frgrd!2 = "lat"
         frgrd!1 = "time"

         frgrd&lon = lon
         frgrd&lat = lat
         frgrd&time = yrs

        end if

        do nt=0, dimsizes(timent)-1
          cyr = datent(nt,0)
          if(any((/1,2,12/).eq.datent(nt,1)))then
            if(datent(nt,1).eq.12.and.ssn.eq."DJF") then
             cyr = cyr+1
            end if
            cssn = "DJF"
          end if
          if(any((/3,4,5/).eq.datent(nt,1)))then
            cssn = "MAM"
          end if
          if(any((/6,7,8/).eq.datent(nt,1)))then
            cssn = "JJA"
          end if
          if(any((/9,10,11/).eq.datent(nt,1)))then
            cssn = "SON"
          end if

          if (cyr.ge.year(0).and.cyr.le.year(1).and.cyr.eq.yr)then
           if (datent(nt,1).ge.mm(0).and.datent(nt,1).le.mm(1))then
            if (datent(nt,2).ge.dd(0).and.datent(nt,2).le.dd(1))then
             if (datent(nt,3).ge.hr(0).and.datent(nt,3).le.hr(1))then

              if (ssn.eq."YYY".or.cssn.eq.ssn)then

                ; for all data points
                ; frgrd(nsel,iyr,:,:) = frgrd(nsel,iyr,:,:)+where(frloc(nt,:,:).gt.0, 1, 0)
                frgrd(nsel,iyr,:,:) = frgrd(nsel,iyr,:,:)+where(frarea(nt,:,:).gt.0, .333, 0)
                ; for NFP
                ; do ifr =1,maxfrN
                ;      Nfr = 0
                ;      Nfr = num(frloc(nt,ind(lat<0),:).eq.ifr)
                ;     ;  print (ifr+"Nfr="+Nfr)
                ;      if (Nfr .gt. 0)then
                ;       ;  print (ifr+ " "+Nfr)
                ;       ;  print (ind(lat.eq.0)+" "+dimsizes(lat))
                ;        do ilat =ind(lat.eq.0),dimsizes(lat)-1
                ;         ;  print (lat(ilat))
                ;          if (num(frloc(nt,ilat,:).eq.ifr).eq.1)then
                ;             ilon =  ind(frloc(nt,ilat,:).eq.ifr)
                ;             ; print(cyr+" "+datent(nt,1)+" "+datent(nt,2)+" "+datent(nt,3)+" "+ifr+" "+lat(ilat)+" "+lon(ilon))
                ;             frgrd(nsel,iyr,ilat,ilon) = frgrd(nsel,iyr,ilat,ilon)+1
                ;             ntstp  = ntstp + 1
                ;             break
                ;          end if
                ;        end do
                ;
                ;      end if
                ;    end do

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
        delete(timent)
        delete(datent)

      end do  ; year
    end do   ;nsel
    ;************************************************
    ; 3. calculate the regression coefficient (slope)
    ;************************************************
    ; rc = new((/dimsizes(lat),dimsizes(lon)/),float) ; add coordinates!
    ;
    ; rc@units = "front/"+ssn
    ;
    ; do ilat =0,dimsizes(lat)-1
    ;   do ilon =0,dimsizes(lon)-1
    ;
    ;     regcoeff = regline( yrs-year(0), frgrd(:,ilat,ilon))          ; slope
    ;     rc(ilat,ilon) = regcoeff
    ;   ;  print(lat(ilat)+" "+lon(ilon)+" "+rc(ilat,ilon) )
    ;   ;  print (rc@yintercept)
    ;
    ;   end do
    ; end do
    ;************************************************
    ; 4. Plotting
    ;************************************************

   print ("Start plotting")
   print ("ntstp="+ntstp)

   ; twist!!!
    ; tmpgrd = frgrd
    ; tmpgrd(:,:299) = frgrd(:,60:)
    ; tmpgrd(:,300:) = frgrd(:,:59)
    ; frgrd = tmpgrd
    ; frgrd&lon = lon
   ;end twist


   ;tstr = cd_calendar(timent(tcyc), 3)

   ;month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
   ;                     "Oct","Nov","Dec"/)
   ;year   = tointeger(datent(tcyc,0))    ; Convert to integer for
   ;month  = tointeger(datent(tcyc,1))    ; use sprinti
   ;day    = tointeger(datent(tcyc,2))
   ;hour   = tointeger(datent(tcyc,3))
   ;date_strm = sprinti("%0.4i", year)+month_abbr(month)+ \
   ;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)
   ;date_str = sprinti("%0.4i", year)+sprinti("%0.2i", month)+ \
  ;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)

   ;print(date_str + "")
   ;status_exit(1)

   ; plot num
   print ("      numbers...")
   avfrgrd = dim_avg_n_Wrap(frgrd,1)

   do nsel =0,1
     if(nsel .eq. 0) then
        yrs = yrs1
     end if
     if(nsel .eq. 1) then
        yrs = yrs2
     end if

      regfrgrd = avfrgrd(nsel,:,:)
      ; regfrgrd = 0
      ; regfrgrd(ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1))) = avfrgrd(nsel,ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1)))
      ; regfrgrd(:,ind(lon.gt.rlon(0).and.lon.le.rlon(1))) = avfrgrd(nsel,:,ind(lon.gt.rlon(0).and.lon.le.rlon(1)))

      wrf_smooth_2d(regfrgrd,3)
      res@tiMainString        = "fronts in "+STR+"  "+sel(nsel)+"  "+ssn+ "  "+ reg

      res@gsnLeftString =str_concat(sprinti("%5.0i", toint(yrs)))
      plot(nsel) = gsn_csm_contour_map_polar(wks,regfrgrd,res)

    end do

  ; draw(plot1)
  ; frame(wks)

  ; plot diff
  print ("      difference..")

  diff = avfrgrd(0,:,:)-avfrgrd(1,:,:)
  wrf_smooth_2d(diff,3)
  copy_VarMeta(regfrgrd,diff)


  resd@tiMainString        = STR+" "+sel(0) +"-"+ sel(1)+", "+ssn+ "  "+ reg

  plot(2) = gsn_csm_contour_map_polar(wks,diff,resd)


  ;************************************************
  ; create panel
  ;************************************************
    resP                    = True                 ; modify the panel plot
    ; resP@gsnPanelMainString = "A common title"     ; new resource added in NCL V6.4.0
    resP@txString           = "fr N: "+STR +" "+sel(0) +"-"+ sel(1)+", "+ ssn + "  "+ regtitle
    gsn_panel(wks,plot,(/2,2/),resP)               ; now draw as one plot

  delete(plot)

  end do ; years
  delete(frgrd)
  delete(yrs)
 end do ; ssn
end do ; reg
end do ; STR
end do ; ihs
end
