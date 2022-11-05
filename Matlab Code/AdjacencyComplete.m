function Yhat = AdjacencyComplete(y,d1,d2,type,radius,sigma,pct,i)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
% Observe 'pct' percent of the entries in Y
strm = RandStream('mt19937ar','Seed',i);
idx = find(rand(strm,d1*d2,1) <= pct/100);
l = length(idx);

%% Define observation model
if type == 1
    % probit/Gaussian noise
    f      = @(x) gausscdf(x,0,sigma);
    fprime = @(x) gausspdf(x,0,sigma);
elseif type == 2
    % Logistic model
    f       = @(x) (1 ./ (1 + exp(-x)));
    fprime  = @(x) (exp(x) ./ (1 + exp(x)).^2);
end

%% Set up optimization problem
options = struct();
options.iterations = 10000;
options.stepMax    = 10000;
options.stepMin    = 1e-4;
options.optTol     = 1e-3;
options.stepMax    = 1e9;

funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);

%% Define constraints
% Use nuclear-norm constraint only
funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);

% Use nuclear-norm plus infinity-norm constraints
% funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);

%% Recover estimate Mhat of M
[Mhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
Mhat = reshape(Mhat,d1,d2);

% [U,S,V] = svd(Mhat);
% Mhat_debias = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Project onto actual rank if known

%% Obtain estimate of the original adjacency matrix
% Define threshold
c = length(find(y(idx)==1))/l;
if c == 0
    c = 1/l;
end
fMhat = f(Mhat);
fMhat = fMhat(:);
fMhat = sort(fMhat,'descend');
thre = fMhat(ceil(length(fMhat)*c));

% Obtain estimate of the upper triangle of Y
Yhat = zeros(d1,d2);
Yhat((f(Mhat) >= thre)) = 1;

% From upper triangle obtain the symmertrix matrix form
% Diag = diag(diag(Yhat));
Yhat = triu(Yhat,1);
Yhat = Yhat + Yhat';
end