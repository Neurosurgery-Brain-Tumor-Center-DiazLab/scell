function [R,rsq]=glm_reg(X,Y,offset)
%function [R,rsq]=glm_reg(X,Y,offset)
%
%IN: X - nxp matrix of n samples by p predictor genes
%    Y - nxm matrix of n samples by m response genes
%    offset - vector of offsets, sample scaling factors
%
%OUT: R - raw residuals from fit
%     rsq - proportion of sum of squares explained by the model

if isempty(gcp('nocreate')), parpool; end
m=size(Y,2);
opt=statset('UseParallel',true);
R=zeros(size(Y));
rsq=zeros(m,1);
warning off;
parfor i=1:m
%    [b,FitInfo]=lassoglm(W,Y(:,i),'poisson','Options',opt,'Offset',offset,'CV',4,'Alpha',1e-4);
%    kpt=find(b(:,FitInfo.IndexMinDeviance));
    mdl=fitglm(X,Y(:,i),'linear','Distribution','poisson','DispersionFlag',true,'Offset',offset);
    R(:,i)=exp(mdl.Residuals.LinearPredictor);
    rsq(i)=mdl.Rsquared.Adjusted;
end
    