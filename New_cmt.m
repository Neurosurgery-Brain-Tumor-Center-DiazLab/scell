function [DataCenter]=New_cmt(Data,p)
%Data should be EntityxFeatures and standardised.
%Calculates ALL centers rather than one at a time
[N,M]=size(Data);
if p==1
    DataCenter=median(Data,1);
    return;
elseif p==2
    DataCenter=mean(Data,1);
    return;
elseif N==1
    DataCenter=Data;
    return;
end
Gradient(1,1:M)=0.001;
OnesIndex(1:N,1)=1;
DataCenter = sum(Data,1)./N;
DistanceToDataCenter=sum(abs(Data - DataCenter(OnesIndex,:)).^p);
NewDataCenter=DataCenter+Gradient;
DistanceToNewDataCenter=sum(abs(Data - NewDataCenter(OnesIndex,:)).^p);
Gradient(1,DistanceToDataCenter < DistanceToNewDataCenter) = Gradient(1,DistanceToDataCenter < DistanceToNewDataCenter).*-1;
while true 
    NewDataCenter = DataCenter + Gradient;
    DistanceToNewDataCenter=sum(abs(Data - NewDataCenter(OnesIndex,:)).^p);  
    Gradient(1,DistanceToNewDataCenter>=DistanceToDataCenter)=Gradient(1,DistanceToNewDataCenter>=DistanceToDataCenter).*0.9;
    DataCenter(1,DistanceToNewDataCenter<DistanceToDataCenter)=NewDataCenter(1,DistanceToNewDataCenter<DistanceToDataCenter);
    DistanceToDataCenter(1,DistanceToNewDataCenter<DistanceToDataCenter)=DistanceToNewDataCenter(1,DistanceToNewDataCenter<DistanceToDataCenter);    
    if all(abs(Gradient)<0.0001), break, end;
end