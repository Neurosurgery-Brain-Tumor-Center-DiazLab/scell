function [nS,sf,d,h,pval]=normalize_samples(S,lbls,show_plt)
%function [nS,sf,d,h,pval]=normalize_samples(S,lbls,show_plt)
%
%IN: S is a nXm matrix of raw fragment counts for n genes and m samples
%    lbls is a cell array of m strings labeling the samples
%    show_plt is a boolean, 1 if you want to display plots
%
%OUT: nS are the variance stabalized and normalized counts
%     sf is a m-vector of scaling factors
%     d is a m-vector of maximal distances to the order-statistics of the
%          geometric mean
%     pval is a vector of p-values testing the hypothesis that the ith
%     cell is enriched over the mean
n=size(S,1); m=size(S,2);
%anscombe transform
nS=2*sqrt(S+3/8);
%trimmed means of order stats
gm=geomean(nS')';
[cs,sidx]=sort(gm);
cs=cumsum(cs);
cs=cs/cs(end);
pn=poissrnd(ones(size(cs))*median(median(S(S>0)))); %poisson noise
pn=pn(sidx);
pn=cumsum(pn);
pn=pn/pn(end);
CS=nS(sidx,:);
%trimmed means of concomitants
CS=cumsum(CS);
for i=1:m, CS(:,i)=CS(:,i)/CS(end,i); end
%compute the max distances to the order-stat's trimmed mean
%compute the point at which they occur and the divergence test p-value
ad=[];ati=[];d=[]; didx=[]; pval=[];zvall=[];cil=[];
for i=1:m
    [td,ti]=max(pn-CS(:,i));
    d=[d;td];
    [atd,ati]=max(abs(pn-CS(:,i)));
    if cs(ati)<CS(ati,i),atd=-atd;end
    ad=[ad;atd];
    didx=[didx;ti];
end
[cx,cxi]=max(pn-cs);
[phat,ci]=binofit(cs(cxi)*n,cxi,.01);
disp(ci)
h=zeros(m,1);pval=zeros(m,1);
for i=1:m
    h(i)=CS(cxi,i)>ci(2);
    pval(i)=binocdf(CS(cxi,i)*n,cxi,phat,'upper');
end
%plot lorenz curves
if show_plt
    figure
    set(gcf,'color','w')
    hold on
    t=linspace(0,1,size(S,1));
    idx1=find(~h);%find(ad>0);
    idx2=find(h);%find(ad<0);
    [~,sidx1]=sort(ad(idx1),'descend');
    [~,sidx2]=sort(ad(idx2),'descend');
    bm1=brewermap(13,'RdBu');
    grad1=colorGradient(bm1(end,:),bm1(9,:),length(idx1));
    grad2=colorGradient(bm1(1,:),bm1(6,:),length(idx2));
    grad1(grad1<0)=0;
    grad2(grad2<0)=0;
    for i=length(idx1):-1:1
        plot(t,CS(:,idx1(sidx1(i))),'color',grad1(i,:))
    end
    %plot(t(1:5:end),cs(1:5:end),'k','LineWidth',2);
    for i=1:length(idx2)
        plot(t,CS(:,idx2(sidx2(i))),'r')
    end
    z=ones(size(cs))*2*sqrt(3/8);
    cz=cumsum(z);cz=cz/cz(end);
    plot(t(1:1000:end),cz(1:1000:end),'ko','LineWidth',1.5,'MarkerSize',8);
    md=mad(S(:));
    rdi=randi([1,length(cs)],[1,floor(.25*length(cs))]);
    z=poissrnd(ones(size(cs))*md);
    [~,dmi]=max(d);
    z(rdi)=poissrnd(nS(sidx(rdi),dmi));
    %cz=cumsum(z);cz=cz/cz(end);
    %plot(t(1:1000:end),cz(1:1000:end),'k+','LineWidth',1.5,'MarkerSize',8);
    %rdi=randi([1,length(cs)],[1,floor(.75*length(cs))]);
    %z=poissrnd(ones(size(cs))*md);
    %[~,dmi]=max(d);
    %z(rdi)=nS(sidx(rdi),dmi);
    %cz=cumsum(z);cz=cz/cz(end);
    %plot(t(1:1000:end),cz(1:1000:end),'k^','LineWidth',1.5,'MarkerSize',8);
    plot(t,cs,'g','LineWidth',3)
    ylabel('% of reads','FontSize',18);
    xlabel('% of genome','FontSize',18);
end
%compute the scaling factors and the scaled samples
md=min(didx);
for i=1:m
    sf(i)=(cs(md)/CS(md,i));
    nS(:,i)=nS(:,i)*sf(i);
end
%visualize the thresholding
if show_plt
    q=quantile(d,[0.1,.25]);%color samples by their quantile in the distribution of divergence tests
    qd10=q(1); qd25=q(2);
    for i=1:m
        if d(i)<=qd25, c{i}='y';else, c{i}='b';end
        if d(i)<=qd10, c{i}='r';end
    end
    x=S(:);x=x(x>0);%raw counts in the top panel
    qt=[.1:.1:1];cbl={};for i=1:length(qt),cbl{i}=num2str(qt(i));end
    q=quantile(x,qt);%identify the number of genes tagged above a given threshold
    N1=zeros(m,length(q)-1);%ie. identify dropouts at different expression thresholds
    for i=1:m
       for j=1:length(q)-1
            N1(i,j)=length(find(S(:,i)>q(j)));
       end
        N1(i,:)=N1(i,:)/sum(N1(i,:));
    end
    N1=N1(:,end:-1:1);
    nSt=zeros(size(nS));for i=1:size(nS,2),tx=S(:,i); nSt(find(tx>0),i)=nS(find(tx>0),i);end
    y=nS(:);%normalized counts on the right panel
    q=quantile(y,qt);
    M1=zeros(m,length(q)-1);
    for i=1:m
        for j=1:length(q)-1
            M1(i,j)=length(find(nSt(:,i)>q(j)));
        end
        M1(i,:)=M1(i,:)/sum(M1(i,:));
    end
    M1=M1(:,end:-1:1);
    %figure
    figure
     bm=brewermap(12,'YlOrRd');
    bm=bm(end:-1:1,:);
    colormap(bm);
    set(gcf,'color','w')
    %plot 1
    subplot(2,1,1)
    bar(N1,'stacked');
    st=sum(N1');
    hold on
    yl=ylim;
    plot(find(pval<=0.05),st(find(pval<=0.05))+.05,'y*')
    plot(find(h),st(find(h))+.05,'r*')
    set(gca,'YLim',[0,1.1],'XLim',[0,m+1])
    title(gca,'raw counts','FontSize',18)
    set(gca,'FontSize',18)
    ylabel({'% genes expressed', 'above quantile'});
    set(gca,'FontSize',18)
    cb=colorbar;
    set(cb,'YDir','reverse','YTickLabel',cbl(end:-1:1));
    ylabel(cb,'Expression quantile','FontSize',18);
    set(cb,'FontSize',18);
    if ~isempty(lbls)
        set(gca,'XtickLabel',lbls);
        rotateXLabels(gca,90);
    end
    %plot 2
    subplot(2,1,2)
    bar(M1,'stacked');
    hold on
    st=sum(M1');
    plot(find(pval<=0.05),st(find(pval<=0.05))+.05,'y*')
    plot(find(h),st(find(h))+.05,'r*')   
    set(gca,'YLim',[0,1.1],'XLim',[0,m+1])
    title(gca,'normalized counts','FontSize',18)
    set(gca,'FontSize',18);
    ylabel({'% genes expressed', 'above quantile'});
    if ~isempty(lbls)
        set(gca,'XtickLabel',lbls);
        rotateXLabels(gca,90);
    end
    cb=colorbar;
    set(cb,'YDir','reverse','YTickLabel',cbl(end:-1:1));
    ylabel(cb,'Expression quantile','FontSize',18);
    set(cb,'FontSize',18);
end