function [pnz,zinf_fdr,iod,iod_fdr]=comp_gene_var_stats(d,wtbar)
%function [pnz,zinf_fdr,iod,iod_fdr]=comp_gene_var_stats(d,wtbar)
%
%IN: d - a single-cell sequencing data object
%    wtbar - a boolean, true if a waitbar should be displayed
%
%OUT: pnz - pnz(i) is the percent of nonzero readcounts (across cells) for
%           gene i
%     zinf_fdr - q values for zinfp, via Storey's method with bootstrapped
%           lambda
%           The p-values test a null hypothesis that the data was drawn
%           from a generalized Poisson distribution with no zero-inflation.
%           c.f. Yang, Z. 2010
%     iod - iod(i) is the index of dispersion (var/mean) for gene i
%     iod_fdr - The underlying p-values test a null hypothesis that readcounts from
%           differnt cells, for the ith gene, were all drawn from poisson
%           distributions with the same mean. c.f. Selby, B. 1965the power function for Selby's test statistic is closed
%               form, and is used to perform this FDR estimates

if isempty(gcp('nocreate')), parpool; end
h=waitbar(0.25,'Performing gene variance analysis...');
[m,n]=size(d.counts);
%index of dispersion analysis:
C=log(d.counts+1);%the dispersion test is performed on a log scale
iod=sum(C')+n^2*(sum((C.^(-1))')).^(-1);iod=iod';%Selby's Wald statistic for the index of dispersion
cr=chi2inv(.99,n-1);%critical point to evaluate non-central chi2, to compute power of W, see Selby '65
                    %we set the rejection region (arbitrarily) at 1%
iod_fdr=1-ncx2cdf(cr,n-1,iod,'upper');%false discovery rate of a test of the null hypothesis of no
                            %differential gene expression between cells, at
                            %the 1% significance level
waitbar(0.5,h,'Performing zero-inflation analysis...');
%zero-inflation analysis:
numnz=sum(d.counts'>0)';%number of cells with non-zero read counts, genewise
num0=n-numnz;
pnz=numnz/n;%percent of cells with non-zero read counts, genewise
%method of Selby, test for zero inflation with poisson null hypothesis
%lmbda=poissfit(d.counts')';%fit poisson distribution, under the null hypothesis of no zero inflation
%p0=poisspdf(0,lmbda);
%zinf_stat=((num0-n*p0).^2)./(n*p0.*(1-p0)-n*mean(C')'.*p0.^2);
%d.zinf_stat=zinf_stat;
%d.zinfp=chi2cdf(zinf_stat,1,'upper');%score test for more zeros than expected under a Poisson model
%
%This approach seems to have lower statistical power, and has been
%abandoned in favor of the (slower) method of Yang, which uses a
%generalized Poisson null hypothesis

%method of Yang et al, test zero inflation with generalized poisson null
lmbda=zeros(m,1);%mean of generalized poisson
options=optimoptions('fminunc','Algorithm','quasi-newton','Display','off','TolFun',1e-3);
parfor i=1:m
    y=d.counts(i,:);yb=mean(y);
    if nnz(y)==0
        zinf_stat(i)=Inf;
        continue;
    end
    [lmbda(i),fval]=fminunc(@(x) abs(sum(y.*(y-1)./(x*(yb-y)+yb*y))-n),mean(y),options);
    zinf_stat(i)=(num0(i)*exp(lmbda(i))-n)^2/(n*exp(lmbda(i))-n-n*lmbda(i)*(lmbda(i)+2)/2);
end
zinfp=chi2cdf(zinf_stat,1,'upper');%score test for more zeros than expected under a Poisson model
zinfp(pnz==1)=1;
zinfp(pnz==0)=0;
[~,q]=mafdr(zinfp);%control for false discovery
zinf_fdr=q';
zinf_fdr(pnz==0)=0;
delete(h)