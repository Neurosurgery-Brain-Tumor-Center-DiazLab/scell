function c=write_cpm(fname,d)
%function out=write_cpm(fname,d)
%
%IN: fname is a string holding the name of the file to write
%    d is a structure holding the data
%
%OUT: c is the number of bytes written to fname

c=0;
f=fopen(fname,'w');
c=c+fprintf(f,'Gene');
for i=1:length(d.slbls),c=c+fprintf(f,'\t%s',d.slbls{i});end
c=c+fprintf(f,'\n');
for i=1:length(d.gsymb)
    c=c+fprintf(f,'%s',d.gsymb{i});
    for j=1:size(d.cpm,2)
        c=c+fprintf(f,'\t%g',d.cpm(i,j));
    end
    c=c+fprintf(f,'\n');
end
