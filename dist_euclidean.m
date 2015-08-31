% Calculates the Euclidean distance between vectors [FAST].
%
% Assume X is an m-by-p matrix representing m points in p-dimensional space and Y is an
% n-by-p matrix representing another set of points in the same space. This function
% compute the m-by-n distance matrix D where D(i,j) is the SQUARED Euclidean distance
% between X(i,:) and Y(j,:).  Running time is O(m*n*p).
%
% If x is a single data point, here is a faster, inline version to use:
%   D = sum( (Y - ones(size(Y,1),1)*x).^2, 2 )';
%
% INPUTS
%   X   - m-by-p matrix of m p dimensional vectors 
%   Y   - n-by-p matrix of n p dimensional vectors 
%
% OUTPUTS
%   D   - m-by-n distance matrix
%
% EXAMPLE
%   X=[randn(100,5)]; Y=randn(40,5)+2;
%   D = dist_euclidean( [X; Y], [X; Y] ); im(D)
%
% DATESTAMP
%   29-Sep-2005  2:00pm
%
% See also DIST_CHISQUARED, DIST_EMD

% Piotr's Image&Video Toolbox      Version 1.03
% Written and maintained by Piotr Dollar    pdollar-at-cs.ucsd.edu 
% Please email me if you find bugs, or have suggestions or questions! 
 
function D = dist_euclidean( X, Y )
    if( ~isa(X,'double') || ~isa(Y,'double'))
        error( 'Inputs must be of type double'); end;
    m = size(X,1); n = size(Y,1);  
    Yt = Y';  
    XX = sum(X.*X,2);        
    YY = sum(Yt.*Yt,1);      
    D = XX(:,ones(1,n)) + YY(ones(1,m),:) - 2*X*Yt;