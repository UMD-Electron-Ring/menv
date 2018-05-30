function [X] = GSmatch2target()
% nonlinear/global search algorithm. Run w/ matlab 2012.

% -- Load otimization settings
load 'runtmp';
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% -- make lattmp variable and load initial opt. parameters into X0
X0 = menv_makekappaarray();

%%
% -- new gs optimization
X = X0; % X is values of optimization variables

sf2 = @(x)svoptfunc(x); % objective
problem = createOptimProblem('fmincon',...
    'x0',X0,...
    'lb',-300*ones(1,length(X0)),...
    'ub',300*ones(1,length(X0)),...
    'objective',sf2);
    
gs = GlobalSearch;
tic
[xg fg flg og] = run(gs,problem);
toc

X = xg; % new strength settings

% -- save error contributions in runtmp and clm.soldata;
f =  optfunc( X ); 
runtmp.f = f;
save 'runtmp' runtmp;
global clm
clm.soldata.f = f;

