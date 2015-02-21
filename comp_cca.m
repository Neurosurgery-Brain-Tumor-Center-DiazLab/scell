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

[A,B,r,U,V,stats] = canoncorr(X,Y);
sidx=find(stats.p<alpha);
A=A(:,sidx);B=B(:,sidx);
U=U(:,sidx);V=V(:,sidx);
r=r(sidx);
%update d with a new factor matrix, the significant canonical factors of gns
rdn=sum(mean(corr(X,V).^2));
if show_plot
    fh=figure;
    set(fh,'color','w');
    ax=gca;
    set(gca,'FontSize',18);
    crs=mean(corr(U,Y));
    [scrs,cidx]=sort(crs,'descend');
    bar(scrs(1:min(length(scrs),10)));%plot the top 10 most correlated genes in the factor
    set(ax,'XTickLabel',d.gsymb(gidx(cidx(1:min(length(cidx),10)))));
    title([fstr ' Tenenhaus redundancy = ' num2str(rdn*100) '%'],'FontSize',18);
    rotateXLabels(ax,90);
end


    
