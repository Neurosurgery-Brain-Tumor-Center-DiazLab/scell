[m,n]=size(d.counts);
C=log(d.counts+1);
iod=sum(C')+n^2*(sum((C.^(-1))')).^(-1);iod=iod';%Selby's Wald statistic for the index of dispersion
d.iod=iod;
cr=chi2inv(.99,n-1);%critical point to evaluate non-central chi2, to compute power of W, see Selby '65
iod_fdr=1-ncx2cdf(cr,n-1,iod,'upper');%false discovery rate of a test of the null hypothesis of no
                            %differential gene expression between cells, at
                            %the 1% significance level
d.iod_fdr=iod_fdr;
numnz=sum(d.counts'>0)';%number of cells with non-zero read counts, genewise
num0=n-numnz;
d.pnz=numnz/n;%percent of cells with non-zero read counts, genewise
lmbda=poissfit(C')';%fit poisson distribution, under the null hypothesis of no zero inflation
p0=poisspdf(0,lmbda);
zinf_stat=((num0-n*p0).^2)./(n*p0.*(1-p0)-n*mean(C')'.*p0.^2);
d.zinf_stat=zinf_stat;
%d.zinf_stat(d.pnz==1)=max(d.zinf_stat);%this is strictly for visualization purposes
d.zinfp=chi2cdf(zinf_stat,1,'upper');%score test for more zeros than expected under a Poisson model
d.zinfp(d.pnz==1)=1;
d.zinfp(d.pnz==0)=0;
[~,q]=mafdr(d.zinfp);%control for false discovery
q(d.pnz==1)=1;
d.zinf_fdr=q;
f=figure;
bet=realmin;
idx_pnz=find(d.zinf_fdr<bet);
idx_pnzc=find(d.zinf_fdr>=bet);
idx_iod=find(d.iod_fdr<bet);
idx_iodc=find(d.iod_fdr>=bet);
gidx={};
idx1=find(d.iod_fdr<bet&d.zinf_fdr>=bet);
idx2=find(d.iod_fdr>=bet|d.zinf_fdr<bet);
gidx(idx1)={'good'};
%gidx(intersect(idx_pnzc,idx_iodc))={'2'};
%gidx(intersect(idx_pnz,idx_iodc))={'3'};
%gidx(intersect(idx_pnzc,idx_iod))={'4'};
gidx(idx2)={'bad'};
gscatter(d.pnz,log(d.iod)/max(log(d.iod)),gidx','br','o',6);
axis square
set(gca,'XLim',[0 1],'Ylim',[0,1]);
