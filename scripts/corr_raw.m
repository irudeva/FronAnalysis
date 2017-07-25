
%%Import frN, SAMssn, STRloc/int


sec = { '0_90S','30_90E.0_90S','90_150E.0_90S','150_210E.0_90S','210_285E.0_90S','300_340E.0_90S'};
ssn = { 'YYY','DJF','MAM','JJA','SON'};

cc_STRloc_frN = double.empty(5,0);
cc_STRint_frN = double.empty(5,0);
for i = 1:5
    
    if i ~=2 
        pc_STRloc_frN_mslp65(i) = partialcorri(STRloc(:,i),frN(:,i),mslp65(:,i));
        cc = corrcoef(STRloc(:,i),frN(:,i));
        cc_STRloc_frN(i) = cc(1,2);
        pc_STRint_frN_mslp65(i) = partialcorri(STRint(:,i),frN(:,i),mslp65(:,i));
        cc = corrcoef(STRint(:,i),frN(:,i));
        cc_STRint_frN(i) = cc(1,2);
   else
        nyrs = size(yrs(:,1));
        pc_STRloc_frN_mslp65(i) = partialcorri(STRloc(2:nyrs,i),frN(2:nyrs,i),mslp65(2:nyrs,i));
        cc = corrcoef(STRloc(2:nyrs,i),frN(2:nyrs,i));
        cc_STRloc_frN(i) = cc(1,2);
        pc_STRint_frN_mslp65(i) = partialcorri(STRint(2:nyrs,i),frN(2:nyrs,i),mslp65(2:nyrs,i));
        cc = corrcoef(STRint(2:nyrs,i),frN(2:nyrs,i));
        cc_STRint_frN(i) = cc(1,2);
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
