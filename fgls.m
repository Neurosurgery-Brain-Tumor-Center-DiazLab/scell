function resid=fgls(X,Y)
%function resid=fgls(X,Y)
%
%IN: Y - response 
%    X - predictor
%OUT:
%    resid - residuals
%feasible generalized least sequares
e=1;
[~,sigma,~]=mvregress(X,Y);
[~,~,resid]=mvregress(X,Y,'algorithm','cwls','covar0',sigma);
resid_old=resid;
while e>0.05
    %ordinary least squares to estimate covariance
    [~,sigma,~]=mvregress(X,Y);
    [~,~,resid]=mvregress(X,Y,'algorithm','cwls','covar0',sigma);
    resid_old=resid;
    e=norm(resid_old-resid)
end


