function d=load_fc_out(fname,bld)
%function d=load_fc_out(fname,bld)
%
%IN:fname - is a string holding the full path of featureCounts output file
%   bld - is a string, 'hum', 'mus' will annotate entrez ids
%OUT:d is a structure holding the data in the file fname
%    d.gsymb - gene symbols
%    d.ent - gene entrez ids
%    d.count - matrix of counts, rows are genes (d.count(i,:) are the
%    expression values for gene d.gsymb{i}), cols are samples
%    d.slbls - cell array of strings, d.slbls{i} is the sample label of
%    d.length - transcript length
%    d.cpm - matrix of transcripts per million, NOT transcript length normalized
if strcmp(bld,'hum')
    load symb2entrez_hsa.mat
    symb2ent=symb2ent_hsa;
elseif strcmp(bld,'mus')
    load symb2entrez_mmu.mat
    symb2ent=symb2ent_mmu;
end
try
    f=fopen(fname);
catch me
    keyboard()
end
dlm=char(9);
t=fgets(f);
t=fgets(f);
nsmp = numel(strfind(t,dlm)) + 1-6;
fclose(f);
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
    if symb2ent.isKey(d.gsymb{i})
        d.ent(i)=symb2ent(d.gsymb{i});
    else
        d.ent(i)=-1;
    end
end
d.length=D{2};
for i=1:length(D)-2
    d.counts(:,i)=D{2+i};
    d.cpm(:,i)=d.counts(:,i)/sum(d.counts(:,i))*1e6;
end


