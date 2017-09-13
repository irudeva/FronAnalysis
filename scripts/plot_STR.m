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
    for iv = 1:2
        ifig = ifig +1;
        if iv == 1
            var = STRloc
        else
            var = STRint
        end

        subplot(nssn,2,ifig)
        for isec = 1:nsec
            isec
            if isec == nsec ||  isec == 2
              p(isec) = plot(squeeze(yrs(1,:)), squeeze(var(isec,issn,:)),':');
            else
              p(isec) = plot(squeeze(yrs(1,:)), squeeze(var(isec,issn,:)),'-');
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
            if iv ==1
                ylim([-42 -24]);
                ylabel('latitude')
            else
                ylim([1013 1027]);
                ylabel('MSLP, hPa')
            end
            if issn == nssn
                xlabel('years')
            end
        end
        if ifig == 10
    %        lgd = legend(flip(charsec))
%            lgd = legend(p,charsecSTR, 'Position',[0.42 0.02 0.2 0.1]);
           columnlegend(6,charsecSTR)
    %        lgd.Location = 'bestoutside'
        end
        %correaltion
%         for isec = 1:nsec
%             i1 = 1;
%             if issn == 2
%                 i1 = 2;
%             end
%             r = corrcoef(STRloc(1,issn,i1:end),STRloc(isec,issn,i1:end));
%             R_STRloc(issn, isec) = r(1,2);
%             stdSTRloc(issn, isec) = std(STRloc(isec,issn,i1:end));
%             clear r
%         end 
%         
    end
    
end
clear charsecSTR var