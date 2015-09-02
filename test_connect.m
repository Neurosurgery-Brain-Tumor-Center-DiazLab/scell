function [c,r]=test_connect(coef,n)
%function [c,r]=test_connect(coef,n)
%
%IN: coef(1,:) x-coordinates of graph
%    coef(2,:) y-coordinates of graph
%    n is the number of "down-samplings" of the edges
%
%OUT: c(i) is the number of connected components after dropping all edges
%>= r(i)
%
%NOTES: for n linearly spaced distances r_i between the maximal and minimal
% edge lengths, edges longer than r_i are dropped and the number of
% connected components is computed.

e=gabrielGraph(coef);
k=max(max(e));
cD=spalloc(k,k,size(e,1));
for i=1:size(e,1)
    cD(e(i,1),e(i,2))=1;
    cD(e(i,2),e(i,1))=1;
end
X=pdist(coef);
r=logspace(log(min(X)),log(max(X)),n);
r=r(end:-1:1);
c(1)=graphconncomp(cD,'DIRECTED',false);
for t=1:length(r)
    for i=1:k-1
        for j=i+1:k
            if X((i-1)*(k-i/2)+j-i)>r(t)
                cD(i,j)=0;
                cD(j,i)=0;
            end
        end
    end 
    c(t+1)=graphconncomp(cD,'DIRECTED',false);
end
r=r(end:-1:1);
r=r/max(r);
c=c/k;
c=c(end:-1:1);