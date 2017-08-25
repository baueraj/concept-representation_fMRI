%Andrew Bauer
%010814

clear all
close all

%% setup

allSubjs = load('./fMRI_data.mat');

% specify which subject(s) to analyze (either individual subj or mean over >1 subj)
allSubjs_analyInd = 1:numel(allSubjs.subjPool);

load ./mammal_trialID_map.mat

%%  

meanSubjRDM_allMasks = nan(30,30,size(allSubjs.mask_pool,1));

for mask_i = 1:size(allSubjs.mask_pool,1)
    mask_ID = char(allSubjs.mask_pool(mask_i, 1));
    
    eval(strcat('allSubjDat = allSubjs.RDMs_allSubj_',mask_ID,';'));
    eval(strcat('allSubj_noFAPassedVox = allSubjs.noFAPassedVox_allSubj_',mask_ID,';'));

    mammalNames = trialIDs_N_names(1:30, 2);
    RDMs_allSubj = allSubjDat(1:30,1:30,allSubjs_analyInd);
    retainSubjInd = ~(isnan(allSubj_noFAPassedVox(allSubjs_analyInd)) | allSubj_noFAPassedVox(allSubjs_analyInd) < 3);
    meanRDM_allSubj = mean(RDMs_allSubj(:,:,retainSubjInd),3);
    
    meanSubjRDM_allMasks(:,:,mask_i) = meanRDM_allSubj;

    %% do clustering over mean of subjects
    clusts = linkage(meanRDM_allSubj, 'average');
    c = cophenet(clusts, meanRDM_allSubj);

    % WITHOUT CUT-OFFS
    fig = figure;    
    [H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right');

    saveas(fig, strcat('./dendrogram_hierarchClust_fMRI_output/',mask_ID,'_S',strrep(num2str(allSubjs_analyInd),' ','')), 'tif');

    %{
    % WITH CUT-OFFS (COLOR THRESHOLDS)    
    fig = figure; 
    incon = inconsistent(clusts);
    T = cluster(clusts, 'cutoff', 0.8);
    k = numel(unique(T));
    temp = sort(clusts(:,3));
    th = temp(size(clusts,1)+2-k);
    [H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right', 'ColorThreshold', th);
        
    saveas(fig, strcat('./dendrogram_hierarchClust_fMRI_output/',mask_ID,'_S',strrep(num2str(allSubjs_analyInd),' ','')), 'tif');
    %}
end

%% do clustering over means of masks
meanRDM_allSubj = mean(meanSubjRDM_allMasks,3);

clusts = linkage(meanRDM_allSubj, 'average');
c = cophenet(clusts, meanRDM_allSubj);

% WITHOUT CUT-OFFS
fig = figure;    
[H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right');

saveas(fig, strcat('./dendrogram_hierarchClust_fMRI_output/meanAllMasks_S',strrep(num2str(allSubjs_analyInd),' ','')), 'tif'); 

%{
% WITH CUT-OFFS (COLOR THRESHOLDS)    
fig = figure; 
incon = inconsistent(clusts);
T = cluster(clusts, 'cutoff', 0.8);
k = numel(unique(T));
temp = sort(clusts(:,3));
th = temp(size(clusts,1)+2-k);
[H, T] = dendrogram(clusts, 0, 'Labels', mammalNames, 'Orientation', 'right', 'ColorThreshold', th);

saveas(fig, strcat('./dendrogram_hierarchClust_fMRI_output/meanAllMasks_S',strrep(num2str(allSubjs_analyInd),' ','')), 'tif');
%}

close all

disp(strcat(mfilename,': done'))