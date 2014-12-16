function [coeff,score,tsq,explained]=visualize_pca(nS,U,dm,ax,pcs)
%function [coeff,score,tsq,explained]=visualize_pca(nS,U,dm,ax,pcs)
%
%IN: nS : matrix of normalized fragment counts, rows==genes, cols==samples
%    U : is a cluster index vector, can be -1 for not
%       clustered
%    dm : is either 2 or 3, representing the dimension in which to plot
%    ax : is the axis to plot into
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
if dm==2
    bm=brewermap(7,'Accent');
    gscatter(coeff(:,1),coeff(:,2),U,bm,'o',10);
    xlabel(['PCA' num2str(pcs(1)) ': ' num2str(explained(1)) '% of variance explained.'],'FontSize',24);
    ylabel(['PCA' num2str(pcs(2)) ': ' num2str(explained(2)) '% of variance explained.'],'FontSize',24);
    title('PCA+MST','FontSize',24);
    axis tight
    set(gcf,'color','w')
else
    cD=squareform(pdist(coeff(:,1:3)));
    [~,ridx]=min(coeff(:,1)); %index of the smallest PCA1 component to use as root
    [T,pred]=graphminspantree(sparse(cD),ridx);
    [un,~,ulocs]=unique(gstr);
    figure
    %bm=brewermap(7,'Accent');
    h=gscatter(coeff(:,1),coeff(:,2),gstr,[],'o',8);
    for i=1:length(h)
        idx=[];
        x=get(h(i),'XData');y=get(h(i),'YData');
        for j=1:length(x)
            idx=[idx;find(coeff(:,1)==x(j)&coeff(:,2)==y(j))];
        end
        set(h(i),'ZData',coeff(idx,3),'MarkerFaceColor',get(h(i),'Color'));
    end
    view(3)
    hold on
    if plt_tr
        for i=1:length(coeff(:,1))
            for j=1:i-1
                if T(i,j)~=0, plot3([coeff(i,1),coeff(j,1)],[coeff(i,2),coeff(j,2)],[coeff(i,3),coeff(j,3)],'c'),end
            end
        end
    end
    xlabel(['PCA1: ' num2str(explained(1)) '% of var.'],'FontSize',24);
    ylabel(['PCA2: ' num2str(explained(2)) '% of var.'],'FontSize',24);
    zlabel(['PCA3: ' num2str(explained(3)) '% of var.'],'FontSize',24);
    title('PCA+MST','FontSize',24);
    axis tight
    set(gcf,'color','w')
end