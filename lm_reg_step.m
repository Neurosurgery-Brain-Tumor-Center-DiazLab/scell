function [var_exp,R]=lm_reg_step(d)
%function [var_exp,R]=lm_reg_step(d)
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

h=waitbar(0,'Normalizing...');
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
%also, X is already on a log scale
for i=1:length(fac_idx)
    X=[X,d.factor{fac_idx(i)}];
    fac_lens(i)=size(d.factor{fac_idx(i)},2);
end
opt=statset('UseParallel',true);
Y=log(d.counts(d.gidx,:)'+1);
var_tot=zeros(size(Y,1),1);
VE=zeros(length(fac_idx),size(Y,2));
R=zeros(size(Y));
k=size(Y,2);
for i=1:k
    waitbar(i/k,h,'Normalizing...');
    mdl=fitlm(X,Y(:,i),'linear');
    R(:,i)=mdl.Residuals.Raw;
    var_tot(i)=mdl.Rsquared.Adjusted;
    vt=zeros(length(fac_idx),1);
    for j=length(fac_idx):-1:2
        terms=['x' num2str(sum(fac_lens(1:(j-1)))+1)];
        for k=(sum(fac_lens(1:(j-1)))+2):sum(fac_lens(1:j))
            terms=[terms '+x' num2str(k)];
        end
        mdl=removeTerms(mdl,terms);
        vt(j)=max(var_tot(i)-mdl.Rsquared.Adjusted,0);
        var_tot(i)=max(var_tot(i)-mdl.Rsquared.Adjusted,0);
    end
    vt(1)=var_tot(i);
    VE(:,i)=vt;
end
delete(h);
for i=1:size(VE,1)
    var_exp{fac_idx(i)}=VE(i,:)';
end


    



    
    
    