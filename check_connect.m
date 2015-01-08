function edst=check_connect(coef,per)
%function edst=check_connect(coef,per)
%
%IN: coef is a nX2 matrix of euclidean coordinates
%    per is a vector of percentages, at which to downsample the ensemble
%
%OUT: edst is a cell array of size length(per), edst{i} is a beta
%distribution object fitted from the edge probabilities down-sampled to a
%percentage of per(i)

edst={};
n=size(coef,1);
for i=1:length(per)
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