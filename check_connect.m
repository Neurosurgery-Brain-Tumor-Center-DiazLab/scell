function [md,h]=check_connect(coef,nit,ep,plt)
%function [md,h]=check_connect(coef,nit,ep,plt)
%
%IN: coef is a nX2 matrix of euclidean coordinates
%    nit is the number of iterations of MCMC to run
%    ep is the amount of perturbation to perform on each sample
%    plt is a bool, if true then will generate a plot and return handle
%d
%OUT: edst is a cell array of size length(per), edst{i} is a beta
%distribution object fitted from the edge probabilities down-sampled to a
%percentage of per(i)
%     h is a plot handle if plt is true, else it is -1

h=-1;
md=[];
n=size(coef,1);
%compute mean for initial input data
e=gabrielGraph(coef);
!rm -f tmp.pairs
f=fopen('tmp.pairs','w');
for j=1:n
    fprintf(f,'%i\t',e(j,1));
    fprintf(f,'%i\n',e(j,2));
end
[status,result]=system('./fitHRG -f tmp.pairs');
x=textscan(result,'%n','Delimiter','\n');
x=x{1}; x=x(end-n:end);
md={mean(x)};
%compute means for perturbed data
parfor i=1:nit
    ptb=(1-2*rand(size(coef)))*ep;
    e=gabrielGraph(coef+ptb);
    tmp=[tempname(pwd) '.pairs'];
    f=fopen(tmp,'w');
    for j=1:n
        fprintf(f,'%i\t',e(j,1));
        fprintf(f,'%i\n',e(j,2));
    end
    [status,result]=system(['./fitHRG -f ' tmp]);
    system(['rm -f ' tmp]);
    x=textscan(result,'%n','Delimiter','\n');
    x=x{1};
    if n>=length(x), md{i+1}=-1;
    else
        x=x(end-n:end);
        md{i+1}=mean(x);
    end
end
for i=1:length(md),y(i)=md{i};end
md=y;
if plt 
    h=boxplot(md)
end
!rm -f tp*