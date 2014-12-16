cD=squareform(pdist(coeff(:,1:2)));
[T,pred]=graphminspantree(sparse(cD),ridx);
hold on
for i=1:length(coeff(:,1))
   for j=1:i-1
       if T(i,j)~=0, plot([coeff(i,1),coeff(j,1)],[coeff(i,2),coeff(j,2)],'c'),end
   end
end