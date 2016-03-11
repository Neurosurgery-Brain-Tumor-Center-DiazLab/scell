function [U,V]=comp_cca(X,Y,alpha)
%function [U,V]=comp_cca(X,Y,alpha)
%
%IN:
%   X and Y are the matrices to correlate, samples-by-genes, Y are the
%   factors that will be regressed out, X are the rest of the genes
%   alpha : significance level for triaging 
%
%OUT:
%   U and V - the matrices of cannonical factors of X and Y resp., which have a p-value <alpha
%       pvalue is computed via Bartlett's approximate chi-squared statistic for H0_K,
%       with Lawley's modification

st=choose_corr_type();
[A,B,r,U,V,stats] = canoncorr2(X,Y,st);
sidx=find(stats.p<alpha);
A=A(:,sidx);B=B(:,sidx);
U=U(:,sidx);V=V(:,sidx);



    
