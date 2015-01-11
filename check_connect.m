function [edst,h]=check_connect(coef,per,plt)
%function [edst,h]=check_connect(coef,per,plt)
%
%IN: coef is a nX2 matrix of euclidean coordinates
%    per is a vector of percentages, at which to downsample the ensemble
%    plt is a bool, if true then will generate a plot and return handle
%
%OUT: edst is a cell array of size length(per), edst{i} is a beta
%distribution object fitted from the edge probabilities down-sampled to a
%percentage of per(i)
%     h is a plot handle if plt is true, else it is -1

h=-1;
edst={};
n=size(coef,1);
m=length(per);
for i=1:m
    idx=randperm(n);
    idx=idx(1:floor(per(i)*n));
    c=coef(sort(idx),:);
    e=gabrielGraph(c);
    !rm -f tmp.pairs
    f=fopen('tmp.pairs','w');
    for j=1:size(e,1)
        fprintf(f,'%i\t',e(j,1));
        fprintf(f,'%i\n',e(j,2));
    end
    [status,result]=system('./fitHRG -f tmp.pairs');
    x=textscan(result,'%n','Delimiter','\n');
    x=x{1}; x=x(end-length(c):end);
    pd=fitdist(x,'beta');
    edst{i}=pd;
end
if plt 
    M=zeros(2,m);
    N=M;
    for i=1:m
        M=paramci(edst{i});
        tb=makedist('beta','a',M(1,1),'b',M(2,2));
        N(1,i)=mean(tb);
        tb=makedist('beta','a',M(2,1),'b',M(1,2));
        N(2,i)=mean(tb);
    end
    h=figure;
    set(h,'color','w');
    bar(N','stacked');
end
    
        
        