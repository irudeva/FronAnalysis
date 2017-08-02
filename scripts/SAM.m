ssn = { 'YYY','DJF','MAM','JJA','SON'};
nssn = length(ssn);

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);

C = {'k','b','r','g','y',[.5 .6 .7],[0 0.75 0.75 ]} 


figure
for issn = 1:nssn
    issn
    subplot(nssn,1,issn)
    for isec = 1:nsec
        isec
        p = plot(squeeze(yrs(1,:)), squeeze(SAMssn(isec,issn,:)),'-');
        p.Color = C{isec};
        if isec == 1
            p.LineWidth = 2
        end 
        hold on;
        title(ssn(issn));
        xlim([yr1-1 yr2+1]);
    end
    if issn == nssn
       lgd = legend(charsec)
       lgd.Location = 'bestoutside'
    end
end
