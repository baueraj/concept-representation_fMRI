%Andrew Bauer
%121714

clear all
close all

%% SETUP: first choose 2 or 3 dimensions, how many cluster for kmeans, and which RDMs dataset
noDim = 3;
noClusters = 8;

allSubjs = load('./fMRI_data.mat'); 

%% specify which subject(s) to analyze (either individual subj or mean over >1 subj)
allSubjs_analyInd = 1:numel(allSubjs.subjPool);

load ./mammal_trialID_map.mat

%%  
 
mammalNames = trialIDs_N_names(1:30,2);
RDMs_allSubj = allSubjs.RDMs_allSubj_nonWLMask_extrasW12_All(1:30,1:30,allSubjs_analyInd);
meanRDM_allSubj = mean(RDMs_allSubj,3);

%% MDS over mean of subjects
options_in = statset('MaxIter',5000,'Display','final');
[scores,stress,disparities] = mdscale(meanRDM_allSubj,noDim,'Start','cmdscale','Replicates',1000,'Options',options_in);
%[scores,stress,disparities] = mdscale(meanRDM_allSubj,noDim,'Start','cmdscale','Replicates',100);

a = mammalNames';
mask_ID = 'nonWL_fromAll';

if noDim == 2

    fig1 = figure;                
    scatter(scores(:,1),scores(:,2));
    dx = 0.005 + rand(1)/150; 
    dy = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,2)+dy, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID); 
    saveas(fig1, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig1, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    [sortX,indX] = sort(scores(:,1),'descend');
    Xranked = [transpose(a(indX)) num2cell(sortX)];

    [sortY,indY] = sort(scores(:,2),'descend');
    Yranked = [transpose(a(indY)) num2cell(sortY)];

    rankedMammalsByDim = {Xranked; Yranked};

elseif noDim == 3

    fig1 = figure;                
    scatter3(scores(:,1),scores(:,2),scores(:,3));
    dx = 0.005 + rand(1)/100;  
    dy = 0.005 + rand(1)/100; 
    dz = 0.005 + rand(1)/100; 
    h = text(scores(:,1)+dx, scores(:,2)+dy, scores(:,3)+dz, a);
    set(h,'Fontsize',14);
    outStub_fig = strcat(mask_ID,'_all3Dims'); 
    saveas(fig1, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');

    fig2 = figure;
    scatter(scores(:,1),scores(:,2));
    dx = 0.005 + rand(1)/150; 
    dy = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,2)+dy, a);
    set(h,'Fontsize',10);        
    outStub_fig = strcat(mask_ID,'_dims-1x2');
    title(strrep(outStub_fig,'_',' '));
    saveas(fig2, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig2, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    fig3 = figure;
    scatter(scores(:,1),scores(:,3));
    dx = 0.005 + rand(1)/150; 
    dz = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,3)+dz, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID,'_dims-1x3'); 
    title(strrep(outStub_fig,'_',' '));
    saveas(fig3, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig3, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    fig4 = figure;
    scatter(scores(:,2),scores(:,3));
    dy = 0.005 + rand(1)/150; 
    dz = 0.005 + rand(1)/150;
    h = text(scores(:,2)+dy, scores(:,3)+dz, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID,'_dims-2x3');
    title(strrep(outStub_fig,'_',' '));
    saveas(fig4, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig4, strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    [sortX,indX] = sort(scores(:,1),'descend');
    Xranked = [transpose(a(indX)) num2cell(sortX)];

    [sortY,indY] = sort(scores(:,2),'descend');
    Yranked = [transpose(a(indY)) num2cell(sortY)];

    [sortZ,indZ] = sort(scores(:,3),'descend');
    Zranked = [transpose(a(indZ)) num2cell(sortZ)];

    rankedMammalsByDim = {Xranked; Yranked; Zranked};
end

% meanRDM_allSubj_tri = triu(meanRDM_allSubj,1);
% eval(strcat('meanRDM_',mask_ID,'=meanRDM_allSubj_tri;'));

outStub_vals = strcat(mask_ID);
save(strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/VALS_',outStub_vals), 'rankedMammalsByDim');

close all

%% kmeans clustering

options_in = statset('MaxIter',5000,'Display','off');
[clustIDVec,clustCentroids,sumd,D] = kmeans(meanRDM_allSubj,noClusters,'distance','correlation','onlinephase','on','emptyaction','drop','replicates',500,'Options',options_in);
silh = silhouette(meanRDM_allSubj, clustIDVec, 'correlation');

clustStore.clustCentroids = clustCentroids;
clustStore.sumd = sumd;
clustStore.D = D;

clustStore.clustIDVec = clustIDVec;
clustStore.silh = silh;

for clust_i = 1:noClusters
    clust_memb_ind = find(clustIDVec == clust_i);
    clust_memb_names = transpose(a(clust_memb_ind));
    eval(strcat('clustStore.clustMembInfo_',num2str(clust_i),' = [num2cell(clust_memb_ind) clust_memb_names];'));
    disp(['clust ' num2str(clust_i)])
    disp('======')
    disp([num2cell(clust_memb_ind) clust_memb_names])
end

disp('======')
disp(['mean silh, ' num2str(noClusters), ' clusters = ' num2str(mean(silh))])
disp('======')
    
save(strcat('./MDS_kMeans_fMRI_output/noDim-',num2str(noDim),'/KMEANS_',outStub_vals), 'clustStore');

disp(strcat(mfilename,': done'))