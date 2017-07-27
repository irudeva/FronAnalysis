
%%Import frN, SAMssn, STRloc/int


sec = { '0_90S','30_90E.0_90S','90_150E.0_90S','150_210E.0_90S','210_285E.0_90S','300_340E.0_90S'};
ssn = { 'YYY','DJF','MAM','JJA','SON'};

yr1= 1979;
yr2 = 2015;
yrs = (yr1:yr2);
nyrs = length(yrs);


cc_STRloc_frN = zeros(6,5);
cc_STRint_frN = zeros(6,5);
for isec = 1:6
for issn = 1:5
    
    if issn ~=2 
        pc_STRloc_frN_mslp65(isec,issn) = partialcorri(squeeze(STRloc(isec,issn,:)),squeeze(frN(isec,issn,:)),squeeze(mslp65(isec,issn,:)));
        cc = corrcoef(STRloc(isec,issn,:),frN(isec,issn,:));
        cc_STRloc_frN(isec,issn) = cc(1,2);
        pc_STRint_frN_mslp65(isec,issn) = partialcorri(squeeze(STRint(isec,issn,:)),squeeze(frN(isec,issn,:)),squeeze(mslp65(isec,issn,:)));
        cc = corrcoef(STRint(isec,issn,:),frN(isec,issn,:));
        cc_STRint_frN(isec,issn) = cc(1,2);
   else
        %nyrs = size(yrs(:,1));
        pc_STRloc_frN_mslp65(isec,issn) = partialcorri(squeeze(STRloc(isec,issn,2:nyrs)),squeeze(frN(isec,issn,2:nyrs)),squeeze(mslp65(isec,issn,2:nyrs)));
        cc = corrcoef(STRloc(isec,issn,2:nyrs),frN(isec,issn,2:nyrs));
        cc_STRloc_frN(isec,issn) = cc(1,2);
        pc_STRint_frN_mslp65(isec,issn) = partialcorri(squeeze(STRint(isec,issn,2:nyrs)),squeeze(frN(isec,issn,2:nyrs)),squeeze(mslp65(isec,issn,2:nyrs)));
        cc = corrcoef(STRint(isec,issn,2:nyrs),frN(isec,issn,2:nyrs));
        cc_STRint_frN(isec,issn) = cc(1,2);
   end
    
end
end

figure
subplot(3,2,1)
plot(1:5,cc_STRloc_frN,'bo');
ylim([0 1]);
xlim([0 6]);
%xlabel('season','FontSize',12,'FontWeight','bold','Color','bl');
xlabel('season','FontSize',12,'FontWeight','bold');
%tickformat('%.2f')
hold on;
plot(1:5,pc_STRloc_frN_mslp65,'bx')
hold on;
plot(1:5,-cc_STRint_frN,'mo');
hold on;
plot(1:5,-pc_STRint_frN_mslp65,'mx')
