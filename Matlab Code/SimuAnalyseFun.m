function [SimuBoots] = SimuAnalyseFun(m,SimuResult,B,pct)
%UNTITLED 此处提供此函数的摘要
% m is Simulation number
% B is Bootstrap sample number
% pct is the pecentage of sampling
SimuBoots = [];
SimuBoots = sparse(SimuBoots);
for k = 1:m
% For each simulation, read the corresponding adjacency matrix
id = find(SimuResult(:,end) == k);
Y = SimuResult(id,1:size(SimuResult,2)-1);

Y(Y == 0) = -1; % change 0 input to -1 input
[d1,d2] = size(Y); % read its dimension

Y_up = triu(Y); % Use the upper triangle of Y only to maintain the symmetricity

y = Y(:); % reshape the matrix into an array

sigma = 0.001; % Noise level

%% Define alpha to be the correct maximum using an oracle
alpha   = 1;
r = rank(Y);
radius  = alpha * sqrt(d1*d2*r);

%% Bootstrap
%B = 1000;

% Define statistics to be Bootstrap-estimated

RelErr =  zeros(1,B); % Relative Error for each Bootstrap sample
Edge = zeros(1,B); % Edge number for each Bootstrap sample
TwoStar = zeros(1,B); % Two-Star number for each Bootstrap sample

% Degree = []; % Degree Matrix to store the degree of each point and its corresponding proportion for each Bootstrap sample

BootsResult = []; % Used to the Bootstrap-estimated Adjacency matrix

for i = 1:B
%% Sample the adjacency matrix
% we need only to obtain the sampled indeces

%pct = 10; % the percentage of matrix to be sampled
%pct = 30; 
%pct = 40;

% Observe 'pct' percent of the entries in Y
strm = RandStream('mt19937ar','Seed',i);
idx = find(rand(strm,d1*d2,1) <= pct/100);
l = length(idx);

%% Define observation model 
% probit/Gaussian noise
f      = @(x) gausscdf(x,0,sigma);
fprime = @(x) gausspdf(x,0,sigma);

% Logistic model
% f       = @(x) (1 ./ (1 + exp(-x)));
% fprime  = @(x) (exp(x) ./ (1 + exp(x)).^2);

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
%funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);

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

Y((Y == -1)) = 0;

%% Compute various statistics
RelErr(i) = norm(Yhat-Y,'fro')/norm(Y,'fro');
Edge(i) = sum(sum(Yhat));
TwoStar(i) = sum(sum(Yhat*Yhat'));

Blabel = i *  ones(d1,1); % Label the result with the index of Bootstrap sample
Yhat = [Yhat Blabel];
BootsResult = vertcat(BootsResult,Yhat);
i
k

end
k
SimuLabel = k * ones(B*d1,1);

BootsResult = [BootsResult SimuLabel];
SimuBoots = vertcat(SimuBoots,BootsResult);
end
end