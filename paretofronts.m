function [M,F] = paretofronts(x,objective,dom,make_plot,parameter)
%PARETOFRONTS from a set of points x with a certain dominance relation dom
%
% [M,F] = PARETOFRONTS(x,objective,dom,make_plot,parameter) can find and/or
% plot the first or all the pareto fronts according to one of 7 relations
% of dominance among points of any dimension. Notice that returning all the
% outputs is considerably more expensive than returning only the first
% output.
%
% Inputs: 
%    x   = table with point in format (N, M) where N is the number of
%    points and M is their dimension
%
%    objective = array specifying if we want to minimize (1) or maximize (0)
%                each dimension (default = 1)
%
%    dom = string or number with a dominace relation (default = pareto
%                                                            dominance)
%          possible values: 1 - 'pareto'
%                           2 - 'lexicographic'
%                           3 - 'extrema'
%                           4 - 'maxdom'
%                           5 - 'cone'
%                           6 - 'epsilon'
%                           7 - 'lorenz'
%
%    make_plot = 0, 1 or 2, plots the points and their fronts
%
%    parameter = available for 'lexicographic' (rank of importance between objectives of length M)
%                              'extrema' (weight vector for each objective of length M)
%                              'epsilon' (epsilon > 0, resolution vector of length M or 1)
%                              'cone' (inclination value lambda (default value is 0.2) of length 1)
%
%
% Outputs:
%    M = list with which elements are in the first front
%    F = list with the front of each element
%
% Examples: 
% 
% members = paretofronts(randn(100,3),[1,1,1],'pareto',1);
%
% [members, fronts] = paretofronts(randn(100,2),[0,0],'epsilon',1,0.5);
%
% members = paretofronts(randn(100,4),[1,1,1,1],'pareto',1);
%
% x = randn(50,2); for i=1:7 subplot(2,4,i); paretofronts(x,[0,0],i,1); end
% 
% See also: randn, rand, plot
%
% $Author: Alan de Freitas $    $Date: 2012/06/07 $    $Revision: 1.0 $
% Copyright: 2012
% http://www.mathworks.com/matlabcentral/fileexchange/authors/255737

if (nargin<1)
    x = [42,58;9,54;60,87;47,26;70,32;70,12;64,94;3,65;7,48;32,64;53,54;65,65;41,54;82,72;72,52;97,99;53,22;33,11;11,11;61,6;78,40;42,45;9,37;27,76;15,63;28,77;44,93;53,97;46,19;88,14;52,70;94,9;64,53;96,53;24,86;68,48;29,39;67,67;70,74;7,52;25,35;22,15;67,59;84,26;34,4;78,75;68,24;1,44;60,69;39,36;92,74;0,39;46,68;42,70;46,44;77,2;32,33;78,42;47,27;4,20;18,82;72,43;47,89;15,39;34,77;61,40;19,81;74,76;24,38;92,22;27,79;77,95;19,33;29,67;9,44;58,83;68,77;55,17;43,86;64,99;65,51;68,88;64,59;95,15;21,20;71,41;24,75;12,83;61,79;45,32;46,53;66,9;77,11;35,14;66,68;42,50;84,19;83,50;26,15;61,5;];
end
if (nargin <2)
    objective = ones(1,size(x,1));
end

% We invert to make 1 for minimization, like standard in other algorithms
% Note that the rest of this code will consider the maximization default
% standard of comparison though
objective = 1 - objective;
   

for i=1:size(x,2)
    if objective(i)==0
        x(:,i) = -x(:,i);
    end
end
if (nargin < 3)
    dom = 1;
end
if ischar(dom)
    switch dom
        case 'pareto'
            dom = 1;
        case 'lexicographic'
            dom = 2;
        case 'extrema'
            dom = 3;
        case 'maxdom'
            dom = 4;
        case 'cone'
            dom = 5;
        case 'epsilon'
            dom = 6;
        case 'lorenz'
            dom = 7;
        otherwise
            disp(['The dominance ',dom,' is not available. Using pareto dominance instead.']);
            dom = 1;
    end
end
if (nargin < 4)
    make_plot = 0;
end
if nargin<5
    if dom == 2
        parameter = 1:size(x,2);
    elseif dom == 3
        parameter = linspace(1,0,size(x,2));
    elseif dom == 5
        parameter = 0.2;
    elseif dom == 6
        parameter = 1;
    else
        parameter = [];
    end
end


n = size(x,1);
if (nargout>1) %If we need all the fronts
    S = zeros(n,n); % Variable for dominance between every 2 elements
end
N = zeros(1,n); % number of elements that dominate p
F = zeros(1,n); % first front
tic;
%% Identifying the first front
%for each p \in P
for p=1:n
    if ((mod(p,50)==1)&&(n>1000))
        if nargout<2
            disp(['Finding first front : ', num2str(p),' of ', num2str(n)]);
        else
            disp(['Finding all fronts : ', num2str(p),' of ', num2str(n)]);
        end
        toc;
        tic;
    end
    %for each q \in P
    for q=randperm(n)
        if dom == 1 % the pareto comparison doesn't call the other functions
        %if (p dominates q) then
            if (sum(x(p,:)>x(q,:))==size(x,2))
                %S(p) := S(p) \union {q} (we put q into S_p, set of solutions dominated by p)
                if (nargout>1) % saves the dominance only if you need many fronts
                    S(p,q) = 1;
                end
            %else, if (p is dominanted by q) then
            elseif (sum(x(q,:)>x(p,:))==size(x,2))
                %n_p = n_p + 1 (increase n_p, or the number of solutions that dominate p)
                N(p) = N(p) + 1;
                if (nargout<2) %breaks if you just need first front
                    break;
                end
            end
        else % the other relations need the functions 
            if (domina(x(p,:),x(q,:),dom,parameter))&&(~(sum(x(q,:)>x(p,:))==size(x,2)))
                %S(p) := S(p) \union {q} (we put q into S_p, set of solutions dominated by p)
                if (nargout>1)
                    S(p,q) = 1;
                end
            %else, if (p is dominanted by q) then
            elseif (domina(x(q,:),x(p,:),dom,parameter))&&(~(sum(x(p,:)>x(q,:))==size(x,2)))
                %n_p = n_p + 1 (increase n_p, or the number of solutions that dominate p)
                N(p) = N(p) + 1;
                if (nargout<2)&&(make_plot==0) %breaks if you just need first front
                    break;
                end
            end
        end
    end
    %if n_p = 0 then (if no element dominates p)
    if (N(p) == 0)
        %F_1 = F_1 \union {p} (p is in the first front)
        F(p) = 1;
    end
    
end

M = F;

%% Other fronts
if (nargout > 1)
    %i=1 begins from the first front
    i=1;
    %while there are elements in front i
    while(sum(F==i)>0)
        %H = empty set
        H = zeros(1,n);
        %for each p \in F_i (for element in the front i)
        for p=1:n
            if (F(p)==i)
                %for each q \in S_p (for each element of S_p, or solutions dominated by p)
                for q=1:n
                    if (S(p,q)==1)
                        %N_q = n_q - 1 (decreases the number of elements dominating q)
                        N(q) = N(q) - 1;
                        %if n_q = 0 (if q is no longer dominated by any point)
                        if (N(q) == 0)
                            %H = H \union {q} (q is in the next front)
                            H(q) = 1;
                        end
                    end
                end
            end
        end
        %i = i + 1, go to the next front
        i = i + 1;
        %F_i = H, the current front is formed by the elements in H
        F = F + H*i;
    end
    %If there is anyone without a front
    if (sum(F==0)>0)
        error('Something wrong happened.');
    end
end

%% Returns the results
%If the user wanted to plot the results
if (make_plot>0)
    if size(x,2) == 2
        for i=1:size(x,2)
            if objective(i)==0
                x(:,i) = -x(:,i);
            end
        end
        hold off;
        if nargout > 1
            for i=max(F):-1:2
                plot(x((F==i),1),x((F==i),2),'.','color',[0 ((-i+max(F))/(max(F)-min(F))) 1-0.4*((-i+max(F))/(max(F)-min(F)))]);
                hold on;
            end
        else
            plot(x((F==0),1),x((F==0),2),'.','color',[0 0 1]);
            hold on;
        end 
        plot(x((M==1),1),x((M==1),2),'.','color',[1 0.3 0]);
        hold on;
        grid on;
        if (make_plot < 2)&&(nargout>1)
            for i=1:n
                text(x(i,1),x(i,2),num2str(F(i)));
            end
        end
        switch dom
            case 1
                title('Pareto dominance');
            case 2
                title(['Lexicographic dominance (Rank =', num2str(parameter),')']);
            case 3
                title(['Extrema dominance (Lambda = ',num2str(parameter),')']);
            case 4
                title('Max dominance');
            case 5
                title(['Cone dominance (Lambda = ',num2str(parameter),')']);
            case 6
                title(['Epsilon dominance (Epsilon = ',num2str(parameter),')']);
            case 7
                title('Lorenz dominance');
            otherwise
                error('Invalid dominance');
        end

    elseif size(x,2) == 3
        for i=1:size(x,2)
            if objective(i)==0
                x(:,i) = -x(:,i);
            end
        end
        hold off; 
        if nargout > 1
        for i=max(F):-1:2
            plot3(x((F==i),1),x((F==i),2),x((F==i),3),'.','color',[0 ((-i+max(F))/(max(F)-min(F))) 1-0.4*((-i+max(F))/(max(F)-min(F)))]);
            hold on;
        end
        else
            plot3(x((F==0),1),x((F==0),2),x((F==0),3),'.','color',[0 0 1]);
            hold on;        
        end
        plot3(x((M==1),1),x((M==1),2),x((M==1),3),'.','color',[1 0.3 0]);
        hold on;
        grid on;
        if (make_plot < 2)&&(nargout>1)
            for i=1:n
                text(x(i,1),x(i,2),x(i,3),num2str(F(i)));
            end
        end
        switch dom
            case 1
                title('Pareto dominance');
            case 2
                title(['Lexicographic dominance (Rank =', num2str(parameter),')']);
            case 3
                title(['Extrema dominance (Lambda = ',num2str(parameter),')']);
            case 4
                title('Max dominance');
            case 5
                title(['Cone dominance (Lambda = ',num2str(parameter),')']);
            case 6
                title(['Epsilon dominance (Epsilon = ',num2str(parameter),')']);
            case 7
                title('Lorenz dominance');
            otherwise
                error('Invalid dominance');
        end
    else
        for i=1:size(x,2)
            if objective(i)==0
                x(:,i) = -x(:,i);
            end
        end
        minimum = min(x);
        for i=1:size(x,1)
            x(i,:) = x(i,:) - minimum;
        end
        maximum = max(x);
        for i=1:size(x,1)
            x(i,:) = x(i,:) ./ maximum;
        end
        hold off;
        plot(x(M==0,:)','--k');
        hold on;
        plot(x(M==1,:)','LineWidth',2);
        switch dom
            case 1
                title('Pareto dominance');
            case 2
                title(['Lexicographic dominance (Rank =', num2str(parameter),')']);
            case 3
                title(['Extrema dominance (Lambda = ',num2str(parameter),')']);
            case 4
                title('Max dominance');
            case 5
                title(['Cone dominance (Lambda = ',num2str(parameter),')']);
            case 6
                title(['Epsilon dominance (Epsilon = ',num2str(parameter),')']);
            case 7
                title('Lorenz dominance');
            otherwise
                error('Invalid dominance');
        end
        set(gca,'ytick',[0,1]);
        set(gca,'xtick',1:size(x,2));
        set(gca,'yticklabel',{'min','max'});
        xlabel('Objectives');
        ylabel('Objective values');
    end
end

end

