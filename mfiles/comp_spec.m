function [odz,stds,sodz]=comp_spec(d,clst,exp_cut,wrt)
%function [odz,stds,sodz]=comp_spec(d,clst,exp_cut,wrt)
%
%IN:
%   d - single-cell data structure
%   clst.slbls - sample ids
%   clst.id - cluster index
%   exp_cut - percentile threshold cutoff to define expressed genes
%   wrt - file name of file to write a TSV file, empty for no output
%
%OUT:
%   odz - length(d.gsymb) -by- length(uniq(clst.id)) log-odds ratios for
%   cluster specificity of each gene
%   std - standard deviations of odz
%   sodz - gene summary specificity score, computed via Der Simion - Laird
%   random effects model

[ctypes,~,uniq_locs]=unique(clst.id);
m=length(d.gsymb);n=length(ctypes);
ncells=size(d.counts,2);
odz=zeros(m,n);%log odds ratio
stds=zeros(m,n);
sodz=zeros(m,1);
for i=1:m
    for j=1:n
        idx=find(uniq_locs==j);
        n11=length(find(d.counts(i,idx)>=exp_cut))+.5;%gene up in ctype{j}
        n10=length(find(d.counts(i,setdiff(1:ncells,idx))>=exp_cut))+.5;%gene up
                                                                 %not in ctype{j}
        n01=length(find(d.counts(i,idx)<exp_cut))+.5;%gene not up in ctype{j}
        n00=length(find(d.counts(i,setdiff(1:ncells,idx))<exp_cut))+.5;%gene
                                                                 %not up,
                                                                 %not in ctype{j}
        odz(i,j)=log((n11*n00)/(n10*n01));
        stds(i,j)=sqrt(sum(1./[n11 n10 n01 n00]));
    end
    w=stds(i,:).^(-2);
    y=odz(i,:);
    ybar=w*y'/sum(w);
    Q=w*((y-ybar).^2)';
    s1=sum(w);s2=sum(w.^2);
    t2=max(0,(Q-n+1)/(s1-s2/s1));
    sodz(i)=sum(y./(stds(i,:).^2+t2))/sum(1./(stds(i,:).^2+t2));
    pnz(i)=sum(d.counts(i,:)>0)/size(d.counts,2);
end
if ~isempty(wrt)
    f=fopen(wrt,'w');
    fprintf(f,'Gene');
    for i=1:n, fprintf(f,'\tCluster_%i',i);end
    fprintf(f,'\tSummary\tPer_nonzero');
    fprintf(f,'\n');
    for i=1:m
        fprintf(f,'%s',d.gsymb{i});
        for j=1:n
            fprintf(f,'\t%g',odz(i,j));
        end
        fprintf(f,'\t%g',sodz(i));
        fprintf(f,'\t%g',pnz(i));
        fprintf(f,'\n');
    end
end