function[y,varargout] = winsor(x,p)
% WINSOR     Winsorize a vector 
% INPUTS   : x - n*1 data vector
%            p - 2*1 vector of cut-off percentiles (left, right) 
% OUTPUTS  : y - winsorized x, n*1 vector
%            i - (optional) n*1 value-replaced-indicator vector
% NOTES    : Let p1 = prctile(x,p(1)), p2 = prctile(x,p(2)). (Note
%            that PRCTILE ignores NaN values). Then 
%            if x(i) < p1, y(i) = min(x(j) | x(j) >= p1)
%            if x(i) > p2, y(i) = max(x(j) | x(j) <= p2)
% EXAMPLE  : x = rand(10,1), y = winsor(x,[10 90])
% AUTHOR   : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 4/15/07


p = prctile(x,p);
i1 = x < p(1); v1 = min(x(~i1));
i2 = x > p(2); v2 = max(x(~i2));
y = x;
y(i1) = v1;
y(i2) = v2;
