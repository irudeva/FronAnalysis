ssn = { 'YYY','DJF','MAM','JJA','SON'};
nssn = length(ssn);

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);

C = {'b','b','r','g','r','k',[.5 .6 .7],[0 0.75 0.75 ]}

charsecSTR = {'\fontname{Times New Roman}\it SH','\it ATL','\it IND','\it AUS','\it WPAC','\it EPAC'}


figure
ifig = 0;
for issn = 1:nssn
    ifig = ifig +1;
    if issn == 2
        ifig = ifig +1;
    end
    
    subplot(3,2,ifig)
    for isec = nsec:-1:1
        isec
        if isec == nsec ||  isec == 2
          p(isec) = plot(squeeze(yrs(1,:)), squeeze(STRloc(isec,issn,:)),':');
        else
          p(isec) = plot(squeeze(yrs(1,:)), squeeze(STRloc(isec,issn,:)),'-');
        end
        -isec+nsec+1
        p(isec).Color = C{-isec+nsec+1};
        p(isec).LineWidth = 1.5
        if isec == 1
            p(isec).LineWidth = 2
        end
        hold on;
        title(ssn(issn));
        xlim([yr1-1 yr2+1]);
        ylim([-42 -24]);
        if any(ismember([1 3 5],ifig))
            ylabel('latitude')
        end
        if any(ismember([5 6],ifig))
            xlabel('years')
        end
    end
    if issn == 1
%        lgd = legend(flip(charsec))
       lgd = legend(p,charsecSTR, 'Position',[0.57 0.7 0.1 0.2]);
%        lgd.Location = 'bestoutside'
    end
    %correaltion
    for isec = 1:nsec
        i1 = 1;
        if issn == 2
            i1 = 2;
        end
        r = corrcoef(STRloc(1,issn,i1:end),STRloc(isec,issn,i1:end));
        R_STRloc(issn, isec) = r(1,2);
        stdSTRloc(issn, isec) = std(STRloc(isec,issn,i1:end));
        clear r
    end 
    
end
clear charsecSTR