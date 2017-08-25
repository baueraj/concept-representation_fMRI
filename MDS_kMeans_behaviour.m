%Andrew Bauer
%072714

clear all
close all

%% SETUP: first choose 2,3, or >3 dimensions, and how many cluster for kmeans
noDim = 3;
noClusters = 7;

%loading the MEAN RDM over subjects in Henley study -- I don't have each
%subject's data
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

%%  

meanOverSubjs = data; %? sqrt(1 - data);
mask_ID = 'origBehav';
  
%% MDS over mean of subjects
options_in = statset('MaxIter',5000,'Display','final');
[scores,stress,disparities] = mdscale(meanOverSubjs,noDim,'Start','cmdscale','Replicates',1000,'Options',options_in);
%[scores,stress,disparities] = mdscale(meanOverSubjs,noDim,'Start','cmdscale','Replicates',100);

a = mammalNames;

if noDim == 2

    fig1 = figure;                
    scatter(scores(:,1),scores(:,2));
    dx = 0.005 + rand(1)/150; 
    dy = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,2)+dy, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID); 
    saveas(fig1, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig1, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    [sortX,indX] = sort(scores(:,1),'descend');
    Xranked = [transpose(a(indX)) num2cell(sortX)];

    [sortY,indY] = sort(scores(:,2),'descend');
    Yranked = [transpose(a(indY)) num2cell(sortY)];

    rankedMammalsByDim = {Xranked; Yranked};

elseif noDim == 3

    fig1 = figure;                
    scatter3(scores(:,1),scores(:,2),scores(:,3));
    dx = 0.001 + rand(1)/100;  
    dy = 0.001 + rand(1)/100; 
    dz = 0.001 + rand(1)/100; 
    h = text(scores(:,1)+dx, scores(:,2)+dy, scores(:,3)+dz, a);
    set(h,'Fontsize',14);
    title('behavior', 'Fontsize', 16);
    outStub_fig = strcat(mask_ID,'_all3Dims'); 
    saveas(fig1, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');

    fig2 = figure;
    scatter(scores(:,1),scores(:,2));
    dx = 0.005 + rand(1)/150; 
    dy = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,2)+dy, a);
    set(h,'Fontsize',10);        
    outStub_fig = strcat(mask_ID,'_dims-1x2');
    title(strrep(outStub_fig,'_',' '));
    saveas(fig2, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig2, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    fig3 = figure;
    scatter(scores(:,1),scores(:,3));
    dx = 0.005 + rand(1)/150; 
    dz = 0.005 + rand(1)/150;
    h = text(scores(:,1)+dx, scores(:,3)+dz, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID,'_dims-1x3'); 
    title(strrep(outStub_fig,'_',' '));
    saveas(fig3, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig3, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    fig4 = figure;
    scatter(scores(:,2),scores(:,3));
    dy = 0.005 + rand(1)/150; 
    dz = 0.005 + rand(1)/150;
    h = text(scores(:,2)+dy, scores(:,3)+dz, a);
    set(h,'Fontsize',10);
    outStub_fig = strcat(mask_ID,'_dims-2x3');
    title(strrep(outStub_fig,'_',' '));
    saveas(fig4, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/FIG_',outStub_fig), 'fig');            
    saveas(fig4, strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/TIF_',outStub_fig), 'tif');

    [sortX,indX] = sort(scores(:,1),'descend');
    Xranked = [transpose(a(indX)) num2cell(sortX)];

    [sortY,indY] = sort(scores(:,2),'descend');
    Yranked = [transpose(a(indY)) num2cell(sortY)];

    [sortZ,indZ] = sort(scores(:,3),'descend');
    Zranked = [transpose(a(indZ)) num2cell(sortZ)];

    rankedMammalsByDim = {Xranked; Yranked; Zranked};
    UNrankedMammalsByDim = [a' num2cell(scores)];

elseif noDim > 3

    rankedMammalsByDim = [];        
    for di = 1:noDim
        [sortD,indD] = sort(scores(:,di),'descend');
        Dranked = [transpose(a(indD)) num2cell(sortD)];

        rankedMammalsByDim = [rankedMammalsByDim; {Dranked}];
    end
    UNrankedMammalsByDim = [a' num2cell(scores)];
end

outStub_vals = strcat(mask_ID);
if noDim > 3
    save(strcat('./MDS_kMeans_behaviour_output/noDim_exceeds3/VALS_',outStub_vals), 'rankedMammalsByDim', 'UNrankedMammalsByDim');
else
    save(strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/VALS_',outStub_vals), 'rankedMammalsByDim', 'UNrankedMammalsByDim');
end

close all

outStub_vals = strcat(mask_ID);

%% kmeans clustering

options_in = statset('MaxIter',5000,'Display','off');
[clustIDVec,clustCentroids,sumd,D] = kmeans(meanOverSubjs,noClusters,'distance','correlation','onlinephase','on','emptyaction','drop','replicates',500,'Options',options_in);
silh = silhouette(meanOverSubjs, clustIDVec, 'correlation');

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

if noDim > 3
    save(strcat('./MDS_kMeans_behaviour_output/noDim_exceeds3/KMEANS_',outStub_vals), 'clustStore');
else
    save(strcat('./MDS_kMeans_behaviour_output/noDim-',num2str(noDim),'/KMEANS_',outStub_vals), 'clustStore');
end

disp(strcat(mfilename,': done'))