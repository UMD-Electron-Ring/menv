function newX = match2period( )
load 'runtmp'

% -- Load otimization settings
maxIter = 40;%round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% -- make lattmp variable
X0 = menv_makekappaarray(); clear X0;

% -- define initial condition
X = [runtmp.ic.x0 runtmp.ic.y0,...
    runtmp.ic.xp0 runtmp.ic.yp0,...
    runtmp.ic.D0 runtmp.ic.Dp0];

% -- least square optimization
options = optimset('LargeScale', 'on', ...
    'Display', 'iter', ...
    'MaxIter', maxIter, ...
    'TolFun', tolFun );
X = lsqnonlin( 'periodfunc', X,[],[],options );

% -- another iterative optimization step to fine-tune
% -- first set all weights equal
[runtmp.weight.xw,runtmp.weight.yw,runtmp.weight.xpw,runtmp.weight.ypw] = deal(1);
save 'runtmp' runtmp

disp(' ');
disp('The X found by lsqnonlin:');
format long;
disp(X);
disp('Continue Iteration ... Residual Error (weight=1)');
for i=1:maxIter
   dx = periodfunc(X);
   X = X+dx/2;
   err = dx*dx';
   disp(X);
   if( err<tolFun )
      break;
   end
end
format;


% -- save error contributions in runtmp and clm.soldata;
f =  periodfunc( X ); 
runtmp.f = f;
save 'runtmp' runtmp;
global clm
clm.soldata.f = f;

% -- return initial conditions
newX = struct();
newX.x0 = X(1);
newX.y0 = X(2);
newX.xp0 = X(3);
newX.yp0 = X(4);
newX.D0 = X(5);
newX.Dp0 = X(6);



