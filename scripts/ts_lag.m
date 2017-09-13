%% Import data from text file.

STR = {'STRlat','STRslp'};
STR_long= {'STR location','(-1)*STR intensity'};

sec = { '0_90S','300_20E.0_90S','30_90E.0_90S','90_150E.0_90S','150_210E.0_90S','210_285E.0_90S'};
%fsec = { '20_40S','300_20E.20_40S','30_90E.20_40S','90_150E.20_40S','150_210E.20_40S','210_285E.20_40S'};
fsec = { '10_90S','300_20E.10_90S','30_90E.10_90S','90_150E.10_90S','150_210E.10_90S','210_285E.10_90S'};
charsec = {'SH','Atl','Ind','Au','WPac','EPac'};
nsec = length(sec);
yr1= 1979;
yr2 = 2015;

% read SAMssn
corrlag = zeros(2,nsec,5,41); 


for iv= 1:2

for isec = 1:nsec

%% Initialize variables.
filename = strcat('/Users/irudeva/work/Projects/Front/FrontAnalysis/output/ts/corrlag.frN_',STR(iv),'.',fsec(isec),'.',int2str(yr1),'_',int2str(yr2),'.txt');
fn = char(filename);
startRow = 2;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%8f%9f%9f%9f%9f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(fn,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
tmp =([dataArray{1:end-1}])';
corrlag(iv,isec,:,:) = tmp(2:end,:);
lag = (tmp(1,:));

end

end
%% Clear temporary variables
clearvars fn filename startRow formatSpec fileID dataArray ans tmp;


ssn = { 'YYY','DJF','MAM','JJA','SON'};
nssn = length(ssn);

C = {'k','b','r','g','y',[.5 .6 .7],[0 0.75 0.75 ]} 



figure
for iv = 1:2
for issn = 1:nssn
    issn
    subplot(nssn,2,(issn-1)*2+iv)
    for isec = 1:nsec
        isec
        if iv==2 
            coeff = -1;
        else
            coeff = 1;
        end
        p(isec) = plot(lag, squeeze(coeff*corrlag(iv,isec,issn,:)),'-');
        grid on;
%         hold on;
%         p2 = plot(lag, squeeze(corrlag(1,isec,issn,:)),'--');
        p(isec).Color = C{isec};
        p(isec).LineWidth = 1.5
%         p2.Color = C{isec};
%         p2.LineWidth = 1.5

        if isec == 1
            p(isec).LineWidth = 2
%             p2.LineWidth = 2
       end 
        hold on;
        if issn ==1
            title(sprintf('%s,   %s', char(STR_long(iv)), char(ssn(issn))))
%             title(strcat(STR_long(iv),',      ',ssn(issn)));
        else
            title(ssn(issn));
        end
        ylim([0 .45]);
        xlim([-130 130]);
    end
    if issn == nssn & iv ==1
       lgd = legend(p(1:3),charsec(1:3));
       lgd.Location = 'northwest'
    end
    if issn == nssn & iv ==2
       lgd = legend(p(4:end),charsec(4:end));
       lgd.Location = 'northeast'
    end
end
end