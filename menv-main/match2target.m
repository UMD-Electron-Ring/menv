function X = match2target( )

% -- this is a stop-gap measure until I figure out what better way to
% manage passing structure variables without hard-disk writes
global runtmp clm
runtmp = clm.tmp.runtmp;


% -- Load otimization settings
%load 'runtmp';
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% -- make lattmp variable and load initial opt. parameters into X0
X0 = menv_makekappaarray();

X = X0;
if( length(X)<=4 ) scale = 'on';
else scale = 'off'; end
options = optimset('LargeScale', scale, ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'optfunc', X,[],[],options );


% -- save error contributions in runtmp and clm.soldata;
f =  optfunc( X ); 
global clm
clm.soldata.f = f;


% -- return optimization variable X

