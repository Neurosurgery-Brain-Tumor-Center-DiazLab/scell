function d=load_fc_out(fname)
%function d=load_fc_out(fname)
%
%IN:fname - is a string holding the full path of featureCounts output file
%   nsmp - the number of samples in the file
%OUT:d is a structure holding the data in the file fname
%    d.gsymb - gene symbols
%    d.entrez - gene entrez ids
%    d.count - matrix of counts, rows are genes (d.count(i,:) are the
%    expression values for gene d.gsymb{i}), cols are samples
%    d.slbls - cell array of strings, d.slbls{i} is the sample label of
%    d.count(:,i)
%    d.length - transcript length
%    d.cpm - matrix of transcripts per million, NOT transcript length normalized
load /Users/aaron/research/fun_genom/data/ent2gsymb_hsa.mat
f=fopen(fname);
dlm=char(9);
t=fgets(f);
t=fgets(f);
nsmp = numel(strfind(t,dlm)) + 1-6;
fclose(f)
keyboard()
f=fopen(fname);
s='';
for i=1:nsmp, s=[s '%s']; end
D=textscan(f,['%*s%*s%*s%*s%*s%*s' s],1,'Headerlines',1);
for i=1:length(D)
    t=D{i}{1};
    idx=strfind(t,'/');
    if ~isempty(idx)
        d.slbls{i}=t(idx(end-1)+1:idx(end)-1);
    else
        d.slbls{i}=t;
    end
end
s='';
for i=1:nsmp, s=[s '%n']; end
D=textscan(f,['%s%*s%*s%*s%*s%n' s]);
d.gsymb=D{1};
for i=1:length(d.gsymb)
    t=min(find(strcmpi(d.gsymb{i},ent2gsymb_hsa.gsymb)));
    if ~isempty(t), d.ent(i)=ent2gsymb_hsa.ent(i); end
end
d.length=D{2};
for i=1:length(D)-2
    d.counts(:,i)=D{2+i};
    d.cpm(:,i)=d.counts(:,i)/sum(d.counts(:,i))*1e6;
end


