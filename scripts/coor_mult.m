%%Import frN, SAMssn, STRloc/int


%sec = { '0_90S','300_340E.0_90S','30_90E.0_90S','90_150E.0_90S','150_210E.0_90S','210_285E.0_90S'};
%charsec = {'SH','Atl','Ind','Au','WPac','EPac'}
ssn = { 'YYY','DJF','MAM','JJA','SON'};

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);


cc_mult = zeros(6,5);
for isec = 1:nsec
for issn = 1:5
    x1 = squeeze(STRloc(isec,issn,:))/10.;
    x2 = -squeeze(STRint(isec,issn,:))/1000.;
    X = [ones(size(x1)) x1 x2 ];
    X1 = [ones(size(x1)) x1 ];
    X2 = [ones(size(x1)) x2 ];
    y = squeeze(frN(isec,issn,:))/1000.;
    [b,bint,r,rint,stats]  = regress(y,X)
    %test
    
    r2_0 = regstats(y,X,'linear','adjrsquare')
    r2_1 = regstats(y,X1,'linear','adjrsquare')
    r2_2 = regstats(y,X2,'linear','adjrsquare')
    %r2_0.adjrsquare
    %r2_1.adjrsquare
    %r2_2.adjrsquare
    if r2_0.adjrsquare < r2_1.adjrsquare || r2_0.adjrsquare < r2_2.adjrsquare
        r2(isec,issn) = 0;
    else
        r2(isec,issn) = 1;
    end
    
    linfit(:,isec,issn) = b(1)+b(2)*x1+b(3)*x2; %+b(4)*x1.*x2
    if issn ~=2 
        c = corrcoef(y,squeeze(linfit(:,isec,issn)));
    else
        c = corrcoef(y(2:nyrs),squeeze(linfit(2:nyrs,isec,issn)));
    end
    cc(isec,issn) = c(1,2);
    
    
    if issn ~=2 
        c = corrcoef(y,x1);
    else
        c = corrcoef(y(2:nyrs),x1(2:nyrs));
    end
    cc_loc(isec,issn) = c(1,2);

    
    if issn ~=2 
        c = corrcoef(y,x2);
    else
        c = corrcoef(y(2:nyrs),x2(2:nyrs));
    end
    cc_int(isec,issn) = c(1,2);

end
end

figure
for issn = 1:nssn
    subplot(nssn,2,(issn-1)*2+1)

     c1 = plot(yrs,squeeze(frN(1,issn,:))/1000. ,'-k');
     c1.LineWidth = 1.5

     hold on;
     c2 = plot(yrs,linfit(:,1,issn) ,'-r');
     c2.LineWidth = 1.5

     %hold on;
     %c3 = plot(yrs,squeeze(STRloc(isec,issn,:)) ,'-b');
     %c3.LineWidth = 1

     %hold on;
     %c4 = plot(yrs,squeeze(STRint(isec,issn,:)) ,'-g');
     %c4.LineWidth = 1
     
     title(ssn(issn));
     xlim([yr1-1 yr2+1]);
     
    if issn == 5
        lgd = legend([c1 c2],{'number of fronts','linear fit'}) ;
        lgd.Location = 'southwest';
        %legend('boxoff')
    end

    subplot(nssn,2,(issn-1)*2+2)
    xleft = .8:5.8;
    xcent = 1:6;
    xright = 1.2:6.2;
    c2 = plot(xleft,cc_loc(:,issn),'b.');
    hold on
    c3 = plot(xright,cc_int(:,issn),'g.');

    c1 = plot(xcent(r2(:,issn)==1),cc(r2(:,issn)==1,issn),'r.','LineWidth',5);
    c4 = plot(xcent(r2(:,issn)==0),cc(r2(:,issn)==0,issn),'rx');
    %title(charsec(isec));
    ylim([0 1]);
    %xlim([0.5 6.5]);
    %grid minor
    %set(gca,'XTick',[0.5:5.5])
    %set(gca,'xgrid','on')
    %set(gca,'XTickLabel',charsec)
    %xlabh = get(gca,'XLabel');
    %set(xlabh,'Position',get(xlabh,'Position') + [10 10 10])

    xlim([0.5 6.5]);
    set(gca,'XTick',[1:6])
    set(gca,'XTickLabel',charsec)
    set(gca,'ticklength',[0 0])
    for iline = 1.5:5.5
        line([iline iline], get(gca, 'ylim'),'Color','black','LineWidth',0.0005);
    end
    if issn == 5
        lgd = legend([c1 c2 c3],{'linear fit','STRloc','-STRint'}) ;
        lgd.Location = 'southeast';
        %legend('boxoff')
    end



end