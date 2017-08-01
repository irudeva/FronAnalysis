
var = {'SH','LH','PE'};

for iv = 1:3
    
    disp(char(var(iv)));
    fnc = char(strcat('/Users/Irina/work/Communication/GhyslaineBoschat/Fluxes/latitudes_when_zm_',var(iv),'_eq_0.nc'));
    vart = ncread(fnc,'time');
    varlat = ncread(fnc,'latitudes');
    
    vartime = (vart)/24 + datenum('1900-01-01 00:00:00');
    clear vart
    
    yr1 = 1979;
    yr2 = 2015;
    
    for i = 1:length(varlat)
       varyr(i,1)=str2double(datestr(vartime(i,1),'yyyy'));
       varmm(i,1)=str2double(datestr(vartime(i,1),'mm')) ;
    end
    
    %YYY
    for iyr = yr1:yr2
        varssn(iyr-yr1+1,1) = mean(varlat(varyr==iyr));
    end
    
    %DJF
    varssn(1,2) = NaN;
    for iyr = yr1+1:yr2
        varssn(iyr-yr1+1,2) = mean(varlat((varyr==iyr-1 & varmm==12) | (varyr==iyr & (varmm==1 | varmm==2))));
    end
    
    
    %MAM
    for iyr = yr1:yr2
        varssn(iyr-yr1+1,3) = mean(varlat(varyr==iyr &(varmm==3 | varmm==4 | varmm==5)));
    end
    
    %JJA
    for iyr = yr1:yr2
        varssn(iyr-yr1+1,4) = mean(varlat(varyr==iyr &(varmm==6 | varmm==7 | varmm==8)));
    end
    
    %SON
    for iyr = yr1:yr2
        varssn(iyr-yr1+1,5) = mean(varlat(varyr==iyr &(varmm==9 | varmm==10 | varmm==11)));
    end
    
    if strcmp(var(iv),'SH') 
        SHssn = varssn;
    else if strcmp(var(iv),'LH') 
        LHssn = varssn;
    else if strcmp(var(iv),'PE') 
        PEssn = varssn;
    end
    end
    end

end

clear fnc i iv iyr var varlat varmm varssn vartime varyr yr1 yr2