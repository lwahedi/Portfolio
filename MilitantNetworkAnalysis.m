%clear
load allianceData.mat
Hatt=Hatt(:,2:75);Hatt(:,67)=[];
rungraph = 0;

%Cluster alliances to create network communities (my level of analysis)
%B=newmangirvan(AlliesMat,439);
%Create Graph
if rungraph ==1
    cols = [1 1 0;1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1];
    colGraph = zeros(size(AlliesMat,1),3);
    counter = 0;
    for i = 1:439; for j = 1:size(B{i},2); colGraph(B{i}(j),:)=cols(counter,:); end; counter = counter+1;if counter==7;counter=1;end;end
    clear counter
    A=drawNetwork(AlliesMat, '-undirected', true, '-nodeColors',colGraph);
end

%ID nonsingleton clusters
clusterInd=zeros(size(B,2),1);
for i = 1:439; if size(B{i},2)>1; clusterInd(i)=1; end; end


%%ID whether cluster is a hub-spoke, and what the density is
hubSpokeClust = zeros(size(B,2),1);
clustDensity = hubSpokeClust;
clustEG = cell(size(hubSpokeClust));
clustEG2= clustEG;
clustEG3= clustEG;
clustEG4= clustEG;
clustSize=hubSpokeClust;
clustLeth=hubSpokeClust;

er=0;
for clust = 1:size(B,2)
    %if a non-singleton
    clustSize(clust) = size(B{clust},2);
    %clustLeth(clust)=sum(leth(B{clust}));
    if size(B{clust},2)>1
        %create subgraph of the cluster, in an adgacency matrix SG
        SG = subgraph(AlliesMat, B{clust});
        %distribution of degrees (temporary, overwritten each round)
        clustDegrees = degrees(SG);
        %The cluster more resembles a hub and spoke the greater the degree
        %of the most connected node relative to the rest. 
        %hubSpokeClust(clust)=(max(zscore(clustDegrees))-median(zscore(clustDegrees)));
        hubSpokeClust(clust)=(max(clustDegrees)-median(clustDegrees))./max(clustDegrees);
        %hubSpokeClust(clust)=mean(zscore(clustDegrees))-median(zscore(clustDegrees));
        %average number of links in entire subgraph
        clustDensity(clust) = link_density(SG);
        [eg, v]=eigencentralityVal(SG);
        if min(eg)<0
            if max(eg<0)
                eg=eg.*-1;
            else
                er=er+1;
            end
        end
        clustEG{clust} = eg.*v;% scaled
        clustEG2{clust} = eg.*clustSize(clust)./861; %by size
        clustEG3{clust}= eg./sqrt(sum(eg.^2)); %euclidean
        clustEG4{clust}= eg; %non-standardized
    end
end

%%Put everything together into dataset
dataMat= zeros(size(AlliesMat,1),4);
egClustAll = zeros(size(AlliesMat,1),4);
%Degree centrality
dataMat(:,1)=degrees(AlliesMat)';

%clusterdensity and hubSpokiness % fixed effects

for clust = 1:size(B,2)
    for g = 1:size(B{clust},2)
        %cluster density
        dataMat(B{clust}(g),2)= clustDensity(clust);
        %hub-spokiness
        dataMat(B{clust}(g),3)=hubSpokeClust(clust);
        %size of cluster
        dataMat(B{clust}(g),4)=clustSize(clust);
        if size(clustEG{clust},1)>0
            %eigenvctor within cluster
            egClustAll(B{clust}(g),1)=clustEG{clust}(g);
            egClustAll(B{clust}(g),2)=clustEG2{clust}(g);
            egClustAll(B{clust}(g),3)=clustEG3{clust}(g);
            egClustAll(B{clust}(g),4)=clustEG4{clust}(g);
        end
    end
end
egCent=eigencentrality(AlliesMat)';

%import other stuff (covariates)
t1=xlsread('Horowitz Potter dataset.xlsx','Sheet1','FN2:FO862');
t2=xlsread('Horowitz Potter dataset.xlsx','Sheet1','FR2:FU862');
t3=xlsread('Horowitz Potter dataset.xlsx','Sheet1','GC2:GC862');
t4=xlsread('Horowitz Potter dataset.xlsx','Sheet1','CM2:CO862');
t1(861,:)=NaN; t2(861,:)=NaN; t3(861,:)=NaN;
dataMat = [dataMat egClustAll(:,2) egClustAll(:,1) t1 t2 t3 t4];
zMat= zscore(dataMat(:,1:6));
%Look only at allied groups (exclude singletons)
nsI=logical(dataMat(:,1)>0);
nsnq=nsI;
nsnq(852)=0;
%skip al qaeda
nQ=ones(861,1);
nQ(852)=0;
nQ=logical(nQ);

%Group level controls
gCont=zeros(size(B,2),6);
clustLeth=zeros(size(clustLeth,1),4); %total, mean, max, sd
nsG=zeros(size(gCont,1),1);%index of non-single clusters
clustCas=clustLeth; %total, mean, max, sd
clustEvents=clustLeth;% total mean max sd
for clust = 1:size(B,2);
    if size(B{clust},2)>1;nsG(clust)=1;end
    gCont(clust,[1 3 4 5 6])=max(dataMat(B{clust},[8 10 11 12 13]));
    gCont(clust,2)=nanmean(dataMat(B{clust},9));
    clustLeth(clust,:)=[sum(dataMat(B{clust},16)),nanmean(dataMat(B{clust},16)),max(dataMat(B{clust},16)),nanstd(dataMat(B{clust},16))];
    clustCas(clust,:)=[sum(sum(dataMat(B{clust},15:16))),nanmean(sum(dataMat(B{clust},15:16),2)),max(sum(dataMat(B{clust},15:16),2)),nanstd(sum(dataMat(B{clust},15:16),2))];
    clustEvents(clust,:)=[sum(dataMat(B{clust},14)),nanmean(dataMat(B{clust},14)),max(dataMat(B{clust},14)),nanstd(dataMat(B{clust},14))];
end
nsG=logical(nsG);

%Create attribute matrix from H&P data
HattA=xlsread('Horowitz Potter dataset.xlsx','Sheet1','BS2:CX862');
HattA=HattA(:,[1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 20 24 28 29 30 31 32]);
HattA(HattA>0)=1;
%%Distance (in tactics employed space)
i=logical(sum(Hatt,2)); i = logical(1-i);
%UseAtt=Hatt;UseAtt(i,:)=NaN;
UseAtt=HattA;
Gdist=cell(size(B));%Average distance from all other members in community
Cdist=zeros(size(B)); %average distance between nodes in community
Cdiverse=zeros(size(B));%distance of community centroid from origin
Gdiverse=0; %distance of group from origin

distMat=zeros(861,2);
i=squareform(pdist([zeros(1,size(UseAtt,2));UseAtt]));
Gdiverse=i(2:size(i,1),1);
distMat(:,2)=Gdiverse;

for clust = 1:size(B,2)
    if size(B{clust})==1
        Gdist{clust}=NaN;
        Cdist(clust)=NaN;
    else
        Cdist(clust)=nanmean(nanmean(squareform(pdist(UseAtt(B{clust},:)))));
        Gdist{clust}=nanmean(squareform(pdist(UseAtt(B{clust},:))),2);
    end
    %Cdiverse, centroid to origin
    i=nanmean(UseAtt(B{clust},:),1); %calculate centroid
    Cdiverse(clust)=pdist([zeros(1,size(UseAtt,2));i]); %distance from zero
    distMat(B{clust},1)=Gdist{clust};
    
end
 
 
% j=zeros(1,size(UseAtt,2));
% for clust = 1:size(B,2)
%     if size(B{clust},2)>1
%         miniAtt=UseAtt(B{clust}',:);
%         dists=squareform(pdist(miniAtt));dists(dists==0)=NaN;
%         Gdist{clust}=nanmean(dists,1);
%         Cdist(clust)=nanmean(nanmean(dists));
%         k1=logical(sum(miniAtt,2)==0);
%         miniAtt(k1,:)=NaN;
%         k=nanmean(miniAtt,1);
%         Cdiverse(clust)=pdist([k;j]);
%         dists=squareform(pdist([j;miniAtt]));
%         Gdiverse{clust}=dists(1,2:size(miniAtt,1)+1);
%     else
%         Gdiverse{clust}=pdist([j;UseAtt(B{clust}',:)]); 
%         if Gdiverse{clust}==0;Gdiverse{clust}=NaN;end;
%         Cdiverse(clust)=pdist([j;UseAtt(B{clust}',:)]);
%         if Gdiverse{clust}==0;Gdiverse{clust}=NaN;end;
%         Gdist{clust}=NaN;
%         Cdist(clust)=NaN;
%     end
% end
% distMat=zeros(861,2);
% for clust = 1:size(B,2)
%     %if size(B,2)>1
%     distMat(B{clust},1)= Gdist{clust}';
%     distMat(B{clust},2)= Gdiverse{clust}';
%     %end
% end


%Time data (for the survival analysis)
YrS=xlsread('Horowitz Potter dataset.xlsx','Sheet1','DB2:DB862');
YrE=xlsread('Horowitz Potter dataset.xlsx','Sheet1','BP2:BP862');
active=xlsread('Horowitz Potter dataset.xlsx','Sheet1','BB2:BB862');
dataMat(:,9)=YrE-YrS;

i=xlsread('Horowitz Potter dataset.xlsx','Sheet1','BG2:BK862');
dataMat(:,10:12)=i(:,[5,2,1]);
clear i


%save('analysisVars.mat', 'dataMat', 'zMat', 'gCont', 'clustLeth',...
%'clustCas','clustEvents', 'egClustAll', 'hubSpokeClust', 'clustDensity',...
%'B', 'nQ','nsnq','nsG')


%%Draw a bunch of figures to look at variable distribution across the
%%network. Preliminary data visualization. 

%Make Eigenvector Graph to show that global EVC is inappropriate 
colGraph(:,1)=zMat(:,6);
colGraph(:,2)=1;
colGraph(:,3)=0;
i=logical(zMat(:,6)>1);
colGraph(i,1)=1;
colGraph(i,2)=0;
colGraph(i,3)=0;
i=logical(zMat(:,6)<0);
colGraph(i,1)=0;
colGraph(i,2)=1-(zMat(i,6).*-3);
colGraph(i,3)=1;
alliesGraph=AlliesMat;
for i=2:861;alliesGraph(1:i-1,i)=0;end;
A=drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);

%Make Eigenvector Graph to show whether local EVC is appropriate 
 colGraph(:,1)=zMat(:,5);
 colGraph(:,2)=1;
 colGraph(:,3)=0;
 i=logical(zMat(:,5)>1);
 colGraph(i,1)=.8;
 colGraph(i,2)=0;
 colGraph(i,3)=0;
 i=logical(zMat(:,5)>2);
 colGraph(i,1)=1;
 i=logical(zMat(:,5)<0);
 colGraph(i,1)=0;
 colGraph(i,2)=1-(zMat(i,5).*-1.5);
 colGraph(i,3)=1;
 alliesGraph=AlliesMat;
 for i=2:861;alliesGraph(1:i-1,i)=0;end;
 A=drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);


%Check whether hubspokes identify hub spoke networks
i=(dataMat(:,3)./5.75);
i=(dataMat(:,3));%./25);
j = logical(i>=.5);
colGraph(j,1)=1;
colGraph(j,2)=1-((i(j)-.5).*2);
colGraph(j,3)=0;
j = logical(i<.5);
colGraph(j,1)=0;
colGraph(j,2)=(i(j).*2);
colGraph(j,3)=1;
A=drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);

%%lethality picture
i = logical(dataMat(:,7)>=1000);
 colGraph(i,1)= 1;
 colGraph(i,2)= 0;
 colGraph(i,3)=0;
i = logical(dataMat(:,7)>=100 & dataMat(:,7)<1000);
 colGraph(i,1)= 1;
 colGraph(i,2)=1-dataMat(i,7)./515;
 i=logical(dataMat(:,7)<100);
 colGraph(i,1)=0;
 colGraph(i,2)=(dataMat(i,7)./60);
 colGraph(i,3)=1;
 i = logical(dataMat(:,7)==0);
 colGraph(i,3)=.7;
 i=isnan(dataMat(:,7));
 colGraph(i,:)=1;
 alliesGraph=AlliesMat;
 for i=2:861;alliesGraph(1:i-1,i)=0;end;
 drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);

%%lethality picture with GTD data
i = logical(leth>=1000);
 colGraph(i,1)= 1;
 colGraph(i,2)= 0;
 colGraph(i,3)=0;
i = logical(leth>=100 & leth<1000);
 colGraph(i,1)= 1;
 colGraph(i,2)=1-leth(i)./1000;
 i=logical(leth<100);
 colGraph(i,1)=0;
 colGraph(i,2)=(leth(i)./100);
 colGraph(i,3)=1;
 i = logical(leth==0);
 colGraph(i,3)=.7;
 i=isnan(leth);
 colGraph(i,:)=1;
 alliesGraph=AlliesMat;
 for i=2:861;alliesGraph(1:i-1,i)=0;end;
 drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);

 
% %  hub-and-spokes
i = logical(dataMat(:,3)>=.75);
 colGraph(i,1)= 1;
 colGraph(i,2)= 0;
 colGraph(i,3)=(1-((dataMat(i,3)-.75).*3));
i = logical(dataMat(:,3)<.75 &dataMat(:,3)>=.5);
 colGraph(i,1)=(dataMat(i,3)-.5).*3;
 colGraph(i,2)=0;
 colGraph(i,3)=1;
i=logical(dataMat(:,3)<.5);
colGraph(i,1:2)=0;
colGraph(i,3)=.7;
 alliesGraph=AlliesMat;
 for i=2:861;alliesGraph(1:i-1,i)=0;end;
 drawNetwork(alliesGraph, '-undirected', true, '-nodeColors',colGraph);
