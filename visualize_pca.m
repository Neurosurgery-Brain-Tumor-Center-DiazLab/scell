function [coeff,score,tsq,explained]=visualize_pca(nS,U,pca_ax,gn_ax,pcs)
%function [coeff,score,tsq,explained]=visualize_pca(nS,U,pca_ax,gn_ax,pcs)
%
%IN: nS : matrix of normalized fragment counts, rows==genes, cols==samples
%    U : is a cluster index vector, can be -1 for not
%       clustered
%    dm : is either 2 or 3, representing the dimension in which to plot
%    pca_ax : is the axis to plot the pca scores into
%    gn_ax : is the axis to plot the gene loadings into
%    pcs : if dm==2, a 2x1 vector of which pc's to plot (e.g. [1,2]) 
%          if dm==3, a 3x1 vector of which pc's to plot
%
%OUT: coeff : PCA coefficeints
%     score : gene loadings
%     tsq : Hotelling's t-square statistic
%     explained : is variance explained

[coeff,score,latent,tsq,explained]=pca(nS);
gstr=cell(length(U),1);
for i=1:length(U)
    if U(i)<0, gstr{i}='Not clustered';
    else, gstr{i}=['Cluster ' num2str(i)];end
end
bm=brewermap(7,'Accent');
%plot the pca scores
axes(pca_ax);
gscatter(coeff(:,1),coeff(:,2),U,bm,'o',10);
xlabel(['PCA' num2str(pcs(1)) ': ' num2str(explained(1)) '% of variance explained.'],'FontSize',20);
ylabel(['PCA' num2str(pcs(2)) ': ' num2str(explained(2)) '% of variance explained.'],'FontSize',20);
title('PCA','FontSize',24);
axis tight;
%plot the gene loadings
axes(gn_ax)
scatter(score(:,1),score(:,2));
xlabel(['PCA' num2str(pcs(1)) ': ' num2str(explained(1)) '% of variance explained.'],'FontSize',20);
ylabel(['PCA' num2str(pcs(2)) ': ' num2str(explained(2)) '% of variance explained.'],'FontSize',20);
title('PCA','FontSize',24);
axis tight;
