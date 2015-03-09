function FinalCentroids = CiK_Means(Data, SmallClusterThreshold, IsDataStandarized, MustLinkIndex, CannotLinkIndex, CannotClusterIndex)
%Constrained iK-means
%For more info see: R. C. de Amorim (2008) Constrained Intelligent K-Means: Improving Results with Limited Previous Knowledge, The 2nd International Conference on Advanced Engineering Computing and Applications in Sciences, IEEE Computer Society Press, pp. 176-180.

%First Step = Standarize data if needed
InitialSize = size(Data,1);
MinEntitiesInCluster = InitialSize * SmallClusterThreshold;
if IsDataStandarized == false
    tmp = Data - repmat(mean(Data), InitialSize ,1);
    Data = tmp./repmat(max(Data) - min(Data),InitialSize , 1);
end
%Second Step = Gets the Sorting Index to be used in the Anomalous Pattern
[desc,SortIndex] = sort(sum(Data.^2,2));

%Third Step = Gets The CloserCannotCluster Array 
CloserCannotCluster = GetCloserCannotCluster(Data, CannotClusterIndex);

%Fourth Step Anomalous Pattern
Centroids = [];
ClusteredEntities = zeros(size(Data,1),1);
TentCentroid = GetCentroid(Data, ClusteredEntities,CannotClusterIndex, SortIndex); % Gets a tentative Centroid     
while ~isempty (TentCentroid)
    while true
        DistDataToCentroid = sum((Data-(repmat(TentCentroid, size(Data,1),1))).^2,2); 
        BelongsToCentroid =  (DistDataToCentroid < sum(Data.^2,2)) &...
            (DistDataToCentroid < CloserCannotCluster);
        %Fix so one entity cannot belong to 2 clusters
        BelongsToCentroid = BelongsToCentroid & ~ClusteredEntities;
        BelongsToCentroid  = SolveCannotLink(BelongsToCentroid, CannotLinkIndex, TentCentroid, Data);
        BelongsToCentroid = SolveMustLink(BelongsToCentroid, MustLinkIndex,ClusteredEntities);
        if size(find(BelongsToCentroid == 1),1) > 1
            NewCentroid = mean(Data(find(BelongsToCentroid==1),:)); %#ok<FNDSB>
        else
            NewCentroid = TentCentroid;
        end
        if isequal(TentCentroid, NewCentroid)
            break
        else
            TentCentroid = NewCentroid;
        end
    end
    Centroids = [Centroids; NewCentroid]; %#ok<AGROW>
    ClusteredEntities(find(BelongsToCentroid==1),1)=size(Centroids,1); %#ok<FNDSB>
    TentCentroid = GetCentroid(Data, ClusteredEntities,CannotClusterIndex,SortIndex); % Gets a tentative Centroid     
end
%remove small centroids
FinalCentroids=[];
for i = 1: size(Centroids,1)
    if size(find(ClusteredEntities==i),1) > MinEntitiesInCluster
        FinalCentroids = [FinalCentroids;Centroids(i,:)]; %#ok<AGROW>
    end
end




function r = GetCentroid(Data, ClusteredEntities,CannotClusterIndex,SortIndex)
%returns a tentative centroid or empty if there isnt one
r = [];
for i = size(Data,1):-1:1
    RealIndex = SortIndex(i,1);
    
    if ClusteredEntities(RealIndex,1) == 0 && isempty(find(CannotClusterIndex==RealIndex,1))
        r = Data(RealIndex,:);
        break;
    end
end


function r = SolveMustLink(BelongsToCentroid, MustLinkIndex,ClusteredEntities)
%The first one to link the other stays
for i = 1 : size(MustLinkIndex,1)
    if BelongsToCentroid(MustLinkIndex(i,1),1) == 1 
        if ClusteredEntities(MustLinkIndex(i,2),1) == 0
            BelongsToCentroid(MustLinkIndex(i,2),1) = 1;
        else
            BelongsToCentroid(MustLinkIndex(i,1),1) = 0;
        end
        
    elseif BelongsToCentroid(MustLinkIndex(i,2),1) == 1 
        if ClusteredEntities(MustLinkIndex(i,1),1) == 0
            BelongsToCentroid(MustLinkIndex(i,1),1) = 1;
        else
            BelongsToCentroid(MustLinkIndex(i,2),1) = 0;
        end
    end
end
r = BelongsToCentroid;



function r = SolveCannotLink(BelongsToCentroid, CannotLinkIndex, TentCentroid, Data)
%If both cannot stay together the closest one stays
for i=1 : size(CannotLinkIndex)
    if BelongsToCentroid(CannotLinkIndex(i,1),1) == 1 && BelongsToCentroid(CannotLinkIndex(i,2),1) == 1
        if sum((Data(CannotLinkIndex(i,1),:) - TentCentroid(1,:)).^2,2) > sum((Data(CannotLinkIndex(i,2),:) - TentCentroid(1,:)).^2,2)
            BelongsToCentroid(CannotLinkIndex(i,1),1) = 0;
        else
            BelongsToCentroid(CannotLinkIndex(i,2),1) = 0;
        end
    end
end
r = BelongsToCentroid;


function CloserCannotCluster = GetCloserCannotCluster(Data, CannotClusterIndex)
%Creates CloserCannotCluster which has the same index as Data and  the
%value is the distance of the closest CannotCluster
CloserCannotCluster = zeros(size(Data,1), 1);
for entity=1:size(Data,1)
    SmallestDistance = sum((Data(entity,:) - Data(CannotClusterIndex(1,1),:)).^2,2);
    for Cannot=2:size(CannotClusterIndex)
        tmpDistance = sum((Data(entity,:) - Data(CannotClusterIndex(Cannot,:),:)).^2,2);
        if tmpDistance < SmallestDistance
            SmallestDistance = tmpDistance;
        end
    end
    CloserCannotCluster(entity,1) = SmallestDistance;
end
   