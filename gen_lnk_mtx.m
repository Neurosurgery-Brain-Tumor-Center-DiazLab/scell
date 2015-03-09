function lnk=gen_lnk_mtx(arch)
%function lnk=gen_lnk_mtx(arch)
%
%IN: arch : a cell array of strings, grouping variable, empty string
%entries have no consraint
%
%OUT: lnk : matrix encoding must- and cant-link rules
%       lnk(i,j)==1 if i and j must link
%       lnk(i,j)==-1 if i and j cant link
%       lnk(i,j)==0 otherwise

[un,~,ulocs]=unique(arch); idx=[]; n=length(un);
for i=2:n, t{i}=find(ulocs==i); idx=union(idx,t{i}); end
lnk=eye(length(arch));
for i=1:n
    lnk(t{i},t{i})=1;
    for j=1:n
        if j~=i
            lnk(t{i},t{j})=-1;
            lnk(t{j},t{i})=-1;
        end
    end
end

