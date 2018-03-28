function X = match2target( )

% -- Load otimization settings
load 'runtmp';
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


f =  optfunc( X ); % run just to get error contributions;
runtmp.f = f;
save 'runtmp' runtmp;

% -- return optimization variable X

