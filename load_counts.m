function d=load_counts(fname,bld,typ)
%function d=load_counts(fname,bld,typ)
%
%IN:fname - is a string holding the full path of featureCounts output file
%   bld - is a string, 'hum', 'mus' will annotate entrez ids
%   typ - is a string, indicating the file type to read in, currently
%         supported types are 'fc' for feature_counts or 'ct' for a TSV count
%         table with sample ids on the first line and gene ids in the first
%         column, an integer table
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
dlm=char(9);%char(9)=='\t'
if strcmp(typ,'fc')
    t=fgets(f);
    t=fgets(f);%the second line has the sample ids
    nsmp = numel(strfind(t,dlm)) + 1-6;%not using the first 6 5 fields, only counts
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
end
if strcmp(typ,'ct')
    t=fgets(f);
    nsmp = numel(strfind(t,dlm));%# cols - (gene id column)
    fclose(f);
    f=fopen(fname);
    s='';
    for i=1:nsmp, s=[s '%s']; end
    D=textscan(f,['%*s' s],1);
    for i=1:length(D), d.slbls{i}=D{i}{1}; end
    s='';
    for i=1:nsmp, s=[s '%n']; end
    D=textscan(f,['%s' s]);
end
d.gsymb=D{1};
for i=1:length(d.gsymb)
    if symb2ent.isKey(d.gsymb{i})
        d.ent(i)=symb2ent(d.gsymb{i});
    else
        d.ent(i)=-1;
    end
end
if strcmp(typ,'fc')
    ofst=2;%count table starts in third column
    d.length=D{2}; 
end
if strcmp(typ,'ct'), ofst=1; end
for i=1:length(D)-ofst
    d.counts(:,i)=D{ofst+i};
    d.cpm(:,i)=d.counts(:,i)/sum(d.counts(:,i))*1e6;
end


