function [var_exp,R]=glm_reg_step(d)
%function [var_exp,R]=glm_reg_step(d)
%
%IN: d - is a single-cell tool data structure
%
%OUT: var_exp - is a cell array of vectors, var_exp{i}(j) is the variance
%               explained by the ith factor, for the jth gene,
%               d.gsymb{j}
%     R - normalized counts
%
%Assumes: d has fields - factors : cell array of factor matrices to include
%                                  in the GLM
%                        factor_ids : a cell array of strings, identifying
%                                     each factor

if isempty(gcp('nocreate')), parpool; end
n=length(d.factor_ids);
var_exp=cell(n,1);
fac_idx=[];
t1=find(strcmp('ERCCs',d.factor_ids));
t2=find(strcmp('Background',d.factor_ids));
t3=find(strcmp('Cyclins',d.factor_ids));
t4=find(strcmp('User_list',d.factor_ids));
if isempty(t1)
    if ~isempty(t2) 
        fac_idx(1)=t2;
    end
else
    fac_idx(1)=t1;
    if ~isempty(t2), fac_idx(2)=t2; end
end
if ~isempty(t3), fac_idx(end+1)=t3; end
if ~isempty(t4), fac_idx(end+t)=t4; end
X=[];fac_lens=zeros(size(fac_idx));%keep track of which columns of X correspond to what factor
%note that this is in the order of fac_idx
for i=1:length(fac_idx)
    X=[X,d.factor{fac_idx(i)}];
    fac_lens(i)=size(d.factor{fac_idx(i)},2);
end
opt=statset('UseParallel',true);
Y=d.counts(d.gidx,:)';
var_tot=zeros(size(Y,1),1);
VE=zeros(length(fac_idx),size(Y,2));
R=zeros(size(Y));
parfor i=1:size(Y,2)
    mdl=fitglm(X,Y(:,i),'linear','Distribution','poisson','DispersionFlag',true,'Offset',log(d.sf));
    R(:,i)=exp(mdl.Residuals.LinearPredictor);
    var_tot(i)=mdl.Rsquared.Adjusted;
    vt=zeros(length(fac_idx),1);
    for j=length(fac_idx):-1:2
        terms=[0,zeros(1,sum(fac_lens(1:j)))];
        terms((sum(fac_lens(1:(j-1)))+2):end)=1;
        mdl=removeTerms(mdl,terms);
        vt(j)=var_tot(i)-mdl.Rsquared.Adjusted;
        var_tot(i)=var_tot(i)-mdl.Rsquared.Adjusted;
    end
    vt(1)=var_tot(i);
    VE(:,i)=vt;
end
for i=1:size(VE,1)
    var_exp{fac_idx(i)}=VE(i,:)';
end


    



    
    
    