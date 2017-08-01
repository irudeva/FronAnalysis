%% Import data from spreadsheet

[~, ~, raw] = xlsread('/Users/irudeva/work/Projects/Front/FrontAnalysis/output/fr_STR_ssn/frN_STRint.1979_2015.xls','corr');

ssn =  raw(4:8,1);
sec =  raw(3,2:7);
for i = 1:6
    frN_STRint(i,:) = raw(4:8,i+1);
end 
clearvars  raw;

[~, ~, raw] = xlsread('/Users/irudeva/work/Projects/Front/FrontAnalysis/output/fr_STR_ssn/frN_STRloc.1979_2015.xls','corr');

ssn =  raw(4:8,1);
sec =  raw(3,2:7);
for i = 1:6
    frN_STRloc(i,:) = raw(4:8,i+1);
end 
clearvars  raw;

[~, ~, raw] = xlsread('/Users/irudeva/work/Projects/Front/FrontAnalysis/output/SAM_STR_ENSO/SAMreg_STRlat.1979_2015.xls','corr');

ssn =  raw(4:8,1);
sec =  raw(3,2:7);
for i = 1:6
    SAM_STRloc(i,:) = raw(4:8,i+1);
end 
clearvars  raw;

[~, ~, raw] = xlsread('/Users/irudeva/work/Projects/Front/FrontAnalysis/output/SAM_STR_ENSO/SAMreg_STRslp.1979_2015.xls','corr');

ssn =  raw(4:8,1);
sec =  raw(3,2:7);
for i = 1:6
    SAM_STRint(i,:) = raw(4:8,i+1);
end 
clearvars  raw;
