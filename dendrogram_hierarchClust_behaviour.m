%Andrew Bauer
%101714

clear all
close all

%%  

load ./behavioural_data.mat

fid = fopen('./mammal_names.txt');
count = 0;   
while 1    
    count = count + 1;    
    tline = fgetl(fid);        
    if ~ischar(tline), break, end    
    textArray_mammalNouns(count).text = tline;
end
fclose(fid);
mammalNames = {textArray_mammalNouns.text};

% meanOverSubjs = sqrt(1 - data);
meanOverSubjs = 1 - data;    
        
%% do clustering over mean of subjects

mask_ID = 'origBehav';
clusts = linkage(meanOverSubjs, 'average');
c = cophenet(clusts, meanOverSubjs);

% WITHOUT CUT-OFFS
fig = figure;    
[H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right');
set(gca,'FontSize',11)
box off
set(gca,'XTick',[])
set(gca,'XColor','w')
% set(gca,'YColor','w')
saveas(fig, strcat('./dendrogram_hierarchClust_behaviour_output/',mask_ID,'_S-all'), 'tif');

%{
% WITH CUT-OFFS (COLOR THRESHOLDS)
fig = figure; 
incon = inconsistent(clusts);
T = cluster(clusts, 'cutoff', 0.8);
k = numel(unique(T));
temp = sort(clusts(:,3));
th = temp(size(clusts,1)+2-k);
[H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right', 'ColorThreshold', th);

saveas(fig, strcat('./dendrogram_hierarchClust_behaviour_output/',mask_ID,'_S-all'), 'tif');
%}

close all

disp(strcat(mfilename,': done'))