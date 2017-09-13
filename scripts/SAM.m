ssn = { 'YYY','DJF','MAM','JJA','SON'};
nssn = length(ssn);

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);

C = {'b','b','r','g','r','k',[.5 .6 .7],[0 0.75 0.75 ]}

charsecSAM = {'\fontname{Times New Roman}\it SH','\it ATL','\it IND','\it AUS','\it WPAC','\it EPAC'}


figure
for issn = 1:nssn
    issn
    subplot(nssn,1,issn)
    for isec = nsec:-1:1
        isec
        if isec == nsec ||  isec == 2
          p(isec) = plot(squeeze(yrs(1,:)), squeeze(SAMssn(isec,issn,:)),':');
        else
          p(isec) = plot(squeeze(yrs(1,:)), squeeze(SAMssn(isec,issn,:)),'-');
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
    end
    if issn == nssn
%        lgd = legend(flip(charsec))
       lgd = legend(p,charsecSAM)
       lgd.Location = 'bestoutside'
    end
    %correaltion
    for isec = 1:nsec
        i1 = 1;
        if issn == 2
            i1 = 2;
        end
        r = corrcoef(SAMssn(1,issn,i1:end),SAMssn(isec,issn,i1:end));
        R(issn, isec) = r(1,2);
        stdSAM(issn, isec) = std(SAMssn(isec,issn,i1:end));
        clear r
    end 
    
end
clear charsecSAM

