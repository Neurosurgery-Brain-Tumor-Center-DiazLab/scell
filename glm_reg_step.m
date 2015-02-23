function var_exp=glm_reg_step(d)
%function var_exp=glm_reg_step(d)
%
%IN: d - is a single-cell tool data structure
%
%OUT: var_exp - is a cell array of vectors, var_exp{i}(j) is the variance
%               explained by the model with the first i factors [d.factor{1}...d.factor{i}] for the jth gene,
%               d.gsymb{j}

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
X=[];fac_lens=zeros(size(fac_idx));
for i=1:length(fac_idx)
    X=[X,d.factor{fac_idx(i)}];
    fac_lens(i)=size(d.factor{fac_idx(i)},2);
end


    
    
    