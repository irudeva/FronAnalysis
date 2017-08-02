
%%Import frN, SAMssn, STRloc/int


%sec = { '0_90S','300_340E.0_90S','30_90E.0_90S','90_150E.0_90S','150_210E.0_90S','210_285E.0_90S'};
%charsec = {'SH','Atl','Ind','Au','WPac','EPac'}
nsec = length(sec);

ssn = { 'YYY','DJF','MAM','JJA','SON'};

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);


cc_STRloc_frN = zeros(6,5);
cc_STRint_frN = zeros(6,5);
for isec = 1:nsec
for issn = 1:5
    
    if issn ~=2 
        pc_STRloc_frN_mslp65(isec,issn) = partialcorri(squeeze(STRloc(isec,issn,:)),squeeze(frN(isec,issn,:)),squeeze(mslp65(isec,issn,:)));
        cc = corrcoef(STRloc(isec,issn,:),frN(isec,issn,:));
        cc_STRloc_frN(isec,issn) = cc(1,2);
        pc_STRint_frN_mslp65(isec,issn) = partialcorri(squeeze(STRint(isec,issn,:)),squeeze(frN(isec,issn,:)),squeeze(mslp65(isec,issn,:)));
        cc = corrcoef(STRint(isec,issn,:),frN(isec,issn,:));
        cc_STRint_frN(isec,issn) = cc(1,2);
        cc = corrcoef(STRloc(isec,issn,:),SAMssn(isec,issn,:));
        cc_STRloc_SAM(isec,issn) = cc(1,2);
        cc = corrcoef(STRint(isec,issn,:),SAMssn(isec,issn,:));
        cc_STRint_SAM(isec,issn) = cc(1,2);
        if isec == 1
            cc = corrcoef(STRloc(isec,issn,:),SHssn(:,issn));
            cc_STRloc_SH(isec,issn) = cc(1,2);
            cc = corrcoef(STRloc(isec,issn,:),LHssn(:,issn));
            cc_STRloc_LH(isec,issn) = cc(1,2);
            cc = corrcoef(STRloc(isec,issn,:),PEssn(:,issn));
            cc_STRloc_PE(isec,issn) = cc(1,2);
            
            cc = corrcoef(STRint(isec,issn,:),SHssn(:,issn));
            cc_STRint_SH(isec,issn) = cc(1,2);
            cc = corrcoef(STRint(isec,issn,:),LHssn(:,issn));
            cc_STRint_LH(isec,issn) = cc(1,2);
            cc = corrcoef(STRint(isec,issn,:),PEssn(:,issn));
            cc_STRint_PE(isec,issn) = cc(1,2);
            
            cc = corrcoef(frN(isec,issn,:),SHssn(:,issn));
            cc_frN_SH(isec,issn) = cc(1,2);
            cc = corrcoef(frN(isec,issn,:),LHssn(:,issn));
            cc_frN_LH(isec,issn) = cc(1,2);
            cc = corrcoef(frN(isec,issn,:),PEssn(:,issn));
            cc_frN_PE(isec,issn) = cc(1,2);
        end
   else
        %nyrs = size(yrs(:,1));
        pc_STRloc_frN_mslp65(isec,issn) = partialcorri(squeeze(STRloc(isec,issn,2:nyrs)),squeeze(frN(isec,issn,2:nyrs)),squeeze(mslp65(isec,issn,2:nyrs)));
        cc = corrcoef(STRloc(isec,issn,2:nyrs),frN(isec,issn,2:nyrs));
        cc_STRloc_frN(isec,issn) = cc(1,2);
        pc_STRint_frN_mslp65(isec,issn) = partialcorri(squeeze(STRint(isec,issn,2:nyrs)),squeeze(frN(isec,issn,2:nyrs)),squeeze(mslp65(isec,issn,2:nyrs)));
        cc = corrcoef(STRint(isec,issn,2:nyrs),frN(isec,issn,2:nyrs));
        cc_STRint_frN(isec,issn) = cc(1,2);
        cc = corrcoef(STRloc(isec,issn,2:nyrs),SAMssn(isec,issn,2:nyrs));
        cc_STRloc_SAM(isec,issn) = cc(1,2);
        cc = corrcoef(STRint(isec,issn,2:nyrs),SAMssn(isec,issn,2:nyrs));
        cc_STRint_SAM(isec,issn) = cc(1,2);
        cc = corrcoef(frN(isec,issn,2:nyrs),SAMssn(isec,issn,2:nyrs));
        cc_frN_SAM(isec,issn) = cc(1,2);
        if isec == 1
            cc = corrcoef(STRloc(isec,issn,2:nyrs),SHssn(2:nyrs,issn));
            cc_STRloc_SH(isec,issn) = cc(1,2);
            cc = corrcoef(STRloc(isec,issn,2:nyrs),LHssn(2:nyrs,issn));
            cc_STRloc_LH(isec,issn) = cc(1,2);
            cc = corrcoef(STRloc(isec,issn,2:nyrs),PEssn(2:nyrs,issn));
            cc_STRloc_PE(isec,issn) = cc(1,2);
            
            cc = corrcoef(STRint(isec,issn,2:nyrs),SHssn(2:nyrs,issn));
            cc_STRint_SH(isec,issn) = cc(1,2);
            cc = corrcoef(STRint(isec,issn,2:nyrs),LHssn(2:nyrs,issn));
            cc_STRint_LH(isec,issn) = cc(1,2);
            cc = corrcoef(STRint(isec,issn,2:nyrs),PEssn(2:nyrs,issn));
            cc_STRint_PE(isec,issn) = cc(1,2);
            
            cc = corrcoef(frN(isec,issn,2:nyrs),SHssn(2:nyrs,issn));
            cc_frN_SH(isec,issn) = cc(1,2);
            cc = corrcoef(frN(isec,issn,2:nyrs),LHssn(2:nyrs,issn));
            cc_frN_LH(isec,issn) = cc(1,2);
            cc = corrcoef(frN(isec,issn,2:nyrs),PEssn(2:nyrs,issn));
            cc_frN_PE(isec,issn) = cc(1,2);
        end
   end
    
end
end

figure
for isec = 1:nsec
    isec
    subplot(3,2,isec)
    xleft = .8:4.8;
    xcent = 1:5;
    xright = 1.2:5.2;
    c1 = plot(xleft,cc_STRloc_frN(isec,:),'ko');
    title(charsec(isec));
    xlim([0.5 5.5]);
%     if isec == 1 
%         ylim([-1 1]);
%         line(get(gca, 'xlim'),[0 0], 'Color','black','LineWidth',0.0005);
%     else
        ylim([0 1]);
 %   end
    %set(gca,'XTick',[1.5:4.5])
    set(gca,'xgrid','off')
    
    set(gca,'XTick',[1:5])
    set(gca,'XTickLabel',ssn)
    set(gca,'ticklength',[0 0])
    for iline = 1.5:4.5
        line([iline iline], get(gca, 'ylim'),'Color','black','LineWidth',0.0005);
    end

    
    hold on;
    p1 = plot(xleft,pc_STRloc_frN_mslp65(isec,:),'kx');
    hold on;
    c2 = plot(xright,-cc_STRint_frN(isec,:),'mo');
    hold on;
    p2 = plot(xright,-pc_STRint_frN_mslp65(isec,:),'mx')
    c3 = plot(xleft,-cc_STRloc_SAM(isec,:),'kd');
    c4 = plot(xright,cc_STRint_SAM(isec,:),'md');
    c5 = plot(xcent,-cc_frN_SAM(isec,:),'b*');
    
    if isec ==1
        lcolor =  [0.8500    0.3250    0.0980];
        lcolor = [0.9290    0.6940    0.1250];
        l1 = plot(xleft-0.1,cc_STRloc_SH,'+','Color',lcolor);
        l2 = plot(xleft-0.1,cc_STRloc_LH,'p','Color',lcolor);
        l3 = plot(xleft-0.1,cc_STRloc_PE,'^','Color',lcolor);

        i1 = plot(xright+0.1,-cc_STRint_SH,'r+');
        i2 = plot(xright+0.1,-cc_STRint_LH,'rp');
        i3 = plot(xright+0.1,-cc_STRint_PE,'r^');

        f1 = plot(xcent,cc_frN_SH,'g+');
        f2 = plot(xcent,cc_frN_LH,'gp');
        f3 = plot(xcent,cc_frN_PE,'g^');

        
        % legend
        % Block 1
        % Axes handle 1 (this is the visible axes)
        ah1 = gca;
        % Legend at axes 1
        lgd = legend(ah1,[l1 l2 l3],{'STRloc / SH','STRloc / LH ','STRloc / PE'}) ;
        %lgd.Location = 'northwest';
        lgd.Position = [0.1 0.62 0.1 0.1];

        % Block 2
        % Axes handle 2 (unvisible, only for place the second legend)
        ah2=axes('position',get(gca,'position'), 'visible','off');
        % Legend at axes 2
        lgd = legend(ah2,[i1 i2 i3],{'-STRint / SH','-STRiny / LH ','-STRint / PE'}) ;
        %lgd.Location = 'north';
        lgd.Position = [0.25 0.62 0.1 0.1];

        % Block 3
        % Axes handle 3 (unvisible, only for place the third legend)
        ah3=axes('position',get(gca,'position'), 'visible','off');
        % Legend at axes 3
        lgd = legend(ah3,[f1 f2 f3],{'frN / SH','frN / LH ','frN / PE'}) ;
        %lgd.Location = 'northeast';
        lgd.Position = [0.4 0.62 0.1 0.1];
    end
    
    
    if isec == 5
        lgd = legend([c1 p1 c3],{' STRloc / frN',' STRloc / frN / mslp65','-STRint / SAM'}) ;
        lgd.Location = 'southwest';
        %legend('boxoff')
    end
    if isec == 6
        ah1 = gca;
        lgd = legend(ah1,[c2 p2 c4],{'-STRint / frN','-STRint / frN / mslp65',' STRint / SAM'}) ;
        lgd.Location = 'southeast';
        %legend('boxoff')

        ah2=axes('position',get(gca,'position'), 'visible','off');
        % Legend at axes 2
        lgd = legend(ah2,[c5],{'-frN / SAM'}) ;
        lgd.Location = 'southwest';

    end
end

clear c1 c2 c3 c4 c5 cc i* lgd p1 p2 x* y* l* f1 f2 f3 a*
