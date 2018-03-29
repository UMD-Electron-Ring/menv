function newX = match2period( )
load 'runtmp'

% -- Load otimization settings
maxIter = 40;%round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% -- make lattmp variable
X0 = menv_makekappaarray(); clear X0;

X = [runtmp.ic.x0 runtmp.ic.y0,...
    runtmp.ic.xp0 runtmp.ic.yp0,...
    runtmp.ic.D0 runtmp.ic.Dp0];
% options = optimset('LargeScale', 'on', ...
%     'Display', 'iter', ...
%     'MaxIter', maxIter, ...
%     'TolFun', tolFun );
% X = lsqnonlin( 'periodfunc', X,[],[],options );

% xw = 1;
% yw = 1;
% xpw = 1;
% ypw = 1;
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
   end;
end;
format;



f =  periodfunc( X ); % run just to get error contributions;
runtmp.f = f;
save 'runtmp' runtmp;

% -- return initial conditions
newX = struct();
newX.x0 = X(1);
newX.y0 = X(2);
newX.xp0 = X(3);
newX.yp0 = X(4);
newX.D0 = X(5);
newX.Dp0 = X(6);



