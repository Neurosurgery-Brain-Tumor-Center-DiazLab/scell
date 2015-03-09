function [U, W, Z, UDistToZ, LoopCount] = SubWkMeans (Data, k, Beta, InitialCentroids, InitialW, p, lnk)
%function [U, W, Z, UDistToZ, LoopCount] = SubWkMeans (Data, k, Beta, InitialCentroids, InitialW, p, lnk)
%
%Parameters:
%Data
%      Dataset, format: Entities x Features
%k
%      Total number of clusters in the dataset
%Beta 
%      Weight Exponent
%Initial Centroids (optional)
%      Initial centroids the algorithm should use, format: k x Features. 
%      Optional, if you don't want to use it just write false.
%InitialW (Optional)
%      Initial set of weights the algorithm should use, format: k x
%      Features. Optional, write false if you don't want to use it.
%p
%      Distance Exponent.
%
%lnk is a size(Data,1) X size(Data,1) matrix. lnk(i,j)==1 iff sample i
%and j must be linked, lnk(i,j)==-1 if sample i and j cant be linked,
%lnk(i,j)==0 otherwise (there is no constraint)
%REMARKS
%      If you want to use the MWK-Means as in the paper, the value of Beta
%      and p should be the same.
%
%Outputs
%
%U
%      Cluster Labels. Clearly they may not directly match the dataset
%      labels (you should use a confusion matrix).
%W
%      Final Weights
%Z
%      Final Centroids
%UDistToZ
%      The distance of each entity to its centroid.
%LoopCount
%      The number of loops the algorithm took to converge. The maximum is
%      hard-coded to 500 (variable MaxLoops)
%
%Note: modified by Aaron Diaz, 3/14, from the original SubWkMeans.m (author Dr. Cordeiro de
%Amorim), to include constraints

%M Lines and N Columns
[M,N] = size(Data);
MaxLoops = 10000;

% Step 1
%Get Random Initial Centroids (VALUES) X Number of K
if exist('InitialCentroids','var') && InitialCentroids(1,1)~=false
    Z = InitialCentroids;
else
    rand('twister', sum(100*clock));
    Z=Data(randperm(M)<=k,:);
end

%Generate Initial set of weights FIXED
if exist('InitialW','var') && InitialW(1,1) ~=false
    W = InitialW;
else
    %W = repmat(1/N, k, N);
    W(1:k,1:N)=1/N;
end

Size_Z=size(Z,1);
temp=zeros(M,Size_Z); OnesIndex(1:M,1)=1;
%compute all pairwise distances to centroids, temp(i,j)==dist from sample
%i to centroid j
for c = 1 : Size_Z
    tmp_Z=Z(c,:);
    tmp_W=W(c,:);
    if p=='f', temp(:,c) = dist(Data, tmp_Z(OnesIndex,:), tmp_W(OnesIndex,:));
    else, temp(:,c) = MinkDist(Data, tmp_Z(OnesIndex,:), p, W(c,:),OnesIndex); end
end

%assign samples to nearest centroids, obeying link rules
U = zeros(M, 1);
OldUDistToZ=U;
for i=1:M
    if U(i)==0 %sample has no cluster
        tidx=find(lnk(i,:)>0);
        if length(tidx)==1
            [~,sidx]=sort(temp(tidx,:));
        else
            [~,sidx]=sort(min(temp(tidx,:)));%sort by smallest to largest
        end                             %distance to a centroid over all
                                     %samples in the must-link group,
                                     %over all centroids
        j=1;%find the first cluster we can enter without violating
                  %a cant-link rule
        while j<=length(sidx)&&violation(U,i,sidx(j),lnk), j=j+1; end
        if j<=length(sidx) 
            tidx=find(lnk(i,:)>0);
            U(tidx)=sidx(j);
            OldUDistToZ(tidx)=temp(tidx,sidx(j));
        else
            disp('more constrained samples than clusters, exiting...')
            keyboard()
            U=[];
            return;
        end
    end
end
OldU = U;


if ~exist('p','var'), p = 'f'; end

LoopCount = 0;

while true
    %Step 2
    %Find Initial U that is minimized for the initials Z and W
    [NewUtp1 UDistToZ]= GetNewU (Data, Z, W.^Beta, p,M,lnk,OldU,OldUDistToZ);
    %If there is no alteration in the labels stop
    if ~any(NewUtp1-U)&&LoopCount~=0, break, end;
    
    %if the labes are equal to the previous-previous lables - stop (cycle)
    if ~any(NewUtp1-OldU)&&LoopCount~=0, break; end;
    
    %Step 3
    OldUDistToZ=UDistToZ;
    OldU = U;
    U = NewUtp1;
    %Get New Centroids
    Ztp1 = GetNewZ(Data, U, k,p, Z);
    %If there is no alteration in the centroids stop
    if Ztp1==Z , break, end;
    
    Z = Ztp1;
    %Step 4
    %Update the Weights
    Wtp1 = GetNewW(Data, U, Z, Beta,p,N);
    %if there is no alteration in the weights, stop
    if norm(Wtp1-W)==0 & LoopCount > 20, break; end;
   W = Wtp1;
   LoopCount = LoopCount + 1;
   if LoopCount>=MaxLoops, break; end;
end




function r=dist(a,b,w,p)
%w here is already w.^Beta
switch (nargin)
    case 2
        %the below is fine, in this case b is a scalar
        r = sum(sum((a - b).^2));
    case 3
        r = sum(sum(((a - b).^2).*w,2),2);
    case 4
        if p =='f'
            r = sum(sum(((a - b).^2).*w));
        end
        
end


function [U UDistToZ]= GetNewU (Data, Z, W, p, M, lnk,OldU,OldUDistToZ)
U=OldU; Size_Z=size(Z,1);
temp=zeros(M,Size_Z); OnesIndex(1:M,1)=1;
%compute all pairwise distances to centroids, temp(i,j)==dist from sample
%i to centroid j
for c = 1 : Size_Z
    tmp_Z=Z(c,:);
    tmp_W=W(c,:);
    if p=='f', temp(:,c) = dist(Data, tmp_Z(OnesIndex,:), tmp_W(OnesIndex,:));
    else, temp(:,c) = MinkDist(Data, tmp_Z(OnesIndex,:), p, W(c,:),OnesIndex); end
end
for i=1:M
        tidx=find(lnk(i,:)>0);
        if length(tidx)==1
            [~,sidx]=sort(temp(tidx,:));
        else
            [~,sidx]=sort(min(temp(tidx,:)));%sort by smallest to largest
        end                             %distance to a centroid over all
                                     %samples in the must-link group,
                                     %over all centroids
        j=1;%find the first cluster we can enter without violating
                  %a cant-link rule
        while j<=length(sidx)&&violation(U,i,sidx(j),lnk), j=j+1; end
        if j<=length(sidx) 
            tidx=find(lnk(i,:)>0);
            U(tidx)=sidx(j);
            UDistToZ(tidx)=temp(tidx,sidx(j));
        else
            disp('something is wrong, I cant enter any cluster, exiting...')
            U=[];
            return;
        end
end




function Z = GetNewZ(Data, U, k,p, OldZ)
Z = OldZ;
if p == 'f'
    for l = 1 : k
        Z(l,:) = mean(Data(U==l,:));
    end
else
    for l = 1 : k
        if sum(U==l)>1
            %if there isnt any entity in the cluster - dont change the Z
            Z(l,:)=New_cmt(Data(U==l,:),p);
        end
    end
end

function W = GetNewW(Data, U, Z, Beta,p,N)
k = size(Z,1);
%N = size(Data,2);
D = zeros(k,N);
%W = zeros(k,N);
W=D;

for l = 1 : k
    for j = 1 : N 
         if p=='f'
            D(l,j) = dist(Data(U==l,j), Z(l,j));  
         else
            D(l,j) = sum(abs(Data(U==l,j)- Z(l,j)).^p);
         end
    end
    
end
D = D + mean(mean(D));


%Calculate the actual Weight
if Beta~=1
    exp = 1/(Beta-1);
    for l = 1 : k
        for j = 1 : N
            tmpD=D(l,j);
            %W(l,j)= 1/sum((repmat(D(l,j),1,N)./D(l,:)).^exp);
            W(l,j)= 1/sum((tmpD(1,ones(N,1))./D(l,:)).^exp);
        end
    end
else
    for l = 1 : k
        [~, MinIndex] = min(D(l,:));
        W(l,MinIndex)=1;
    end
end 

function r= MinkDist(x, y, p, w,OnesIndex)
%calculates the  Minkowski distance in between x and y
%r = sum((abs(x - y).^p).*repmat(w,size(x,1),1),2).^(1/p);
%r = sum((abs(x - y).^p).*w(ones(Size_x,1),:),2).^(1/p);
r = sum((abs(x - y).^p).*w(OnesIndex,:),2).^(1/p);

function t=violation(U,idx,c,lnk)
%U is the vector of current cluster assignments, U(i)==0 if sample i has
%not been assigned a cluster, U(i)==c>0 if U(i) is in cluster c
%idx is an index into U, of the sample that we would like to move
%c is the ordinal of the cluster we would like to move into
%lnk is a matrix of must link and cant link rules
%
%This is function checks to see if moving sample idx into cluster c (and
%dragging its must link partners with it) would violate a cant link rule
t=0;
tidx=find(lnk(idx,:)>0);%find sample idx's must link partners
cli=[];
for i=1:length(tidx) %over each must-link partner aggregate the cant link indices
    cli=union(cli,find(lnk(tidx(i),:)<0));
end
t=~isempty(intersect(cli,find(U==c)));
 



