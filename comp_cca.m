function d=comp_cca(d,gns,fstr,alpha,show_plot)
%function d=comp_cca(d,gns,fstr,alpha,show_plot)
%
%IN:
%d : norm_tool's single-cell transcriptomics data structure
%gns : a list of genes from d.gsymb, forming one of the CCA groups
%fstr : an arbitrary string identifying the cannonical factors of gns
%alpha : significance level for triaging 
%show_plot : produce plots
%
%OUT:
%d : d will have fields d.facs and d.facm updated, to include new factors
%    those fields will be added if not there. d.facs is a string array of
%    factor identifiers and d.facm is a cell array of matrices of factors
%    to optionally regress out

n=size(d.counts,1);%number of genes, overall
gidx=[];
for i=1:length(gns)
    t=min(find(strcmp(gns{i},d.gsymb)));
    if ~isempty(t),gidx=[gidx;t];end
end
if isempty(gidx),disp('error no genes found');return;end
X=d.counts(setdiff(1:n,gidx),:)';
Y=d.counts(gidx,:)';
[A,B,r,U,V,stats] = canoncorr(X,Y);
sidx=find(stats.p<alpha);
A=A(:,sidx);B=B(:,sidx);
U=U(:,sidx);V=V(:,sidx);
r=r(sidx);
%update d with a new factor matrix, the significant canonical factors of gns
if ~isfield(d,'facs')
    d.facs={};d.facm={};t=1;
else
    t=min(find(strcmp(fstr,d.facs)));
    if isempty(t),t=length(d.facs)+1;end
end
d.facs{t}=fstr;d.facm{t}=V;
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
%write a ranking of genes in the factor, by mean canonical cross
%correlation
[~,cidx]=sort(abs(crs),'descend');
[fname pname]=uiputfile([fstr '_correlation_rank.tsv'],['Where should I save ' fstr ', ranked by correlation?']);
f=fopen(fullfile(pname,fname),'w');
fprintf(f,[fstr '_gene\tMean_correlation\n']);
for i=1:length(sidx)
    fprintf(f,'%s\t',d.gsymb{gidx(cidx(i))});
    fprintf(f,'%g\n',crs(cidx(i)));
end
fclose(f);
%write a list of top correlated genes with the factor
[fname pname]=uiputfile([fstr '_correlated_genes.tsv'],['Where should I save genes most correlated with' fstr]);
f=fopen(fullfile(pname,fname),'w');
fprintf(f,['Gene\tMean_correlation_with_' fstr 'COV\n']);
crs=max(corr(X,V)');
cut=quantile(crs,.5);
cvs=var(X)./(mean(X).^2);%coefficient of variation
[scvs,cvidx]=sort(cvs,'descend');
cgidx=setdiff(1:n,gidx);
for i=1:length(cvidx)
    fprintf(f,'%s\t',d.gsymb{cgidx(cvidx(i))});
    fprintf(f,'%g\t',crs(cvidx(i)));
    fprintf(f,'%g\n',scvs(i));
end
fclose(f);

    
