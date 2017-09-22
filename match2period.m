function X = match2period( paramfile )
global kx ky K Ex Ey ds;
global KX KY nsteps;
global xw yw xpw ypw;
load( paramfile );

% Global beam parateters
K = runtmp.perveance; % Pervence
Ex = runtmp.emitance; % Emmitance x
Ey = Ex;               % Emmitance y
% Initail conditions
x0 = runtmp.x0;
y0 = runtmp.y0;
xp0 = runtmp.xp0;
yp0 = runtmp.yp0;
% weights
xw = runtmp.xw;
yw = runtmp.yw;
xpw = runtmp.xpw;
ypw = runtmp.ypw;
% Numerical parameters
max_d = runtmp.distance;
min_d = 0.0;
ds = runtmp.stepsize;          % Step-size
nsteps = round((max_d-min_d)/ds)+1;  % steps
% Lattice
ele = runtmp.ele;      % element: 'Q'/'S'
loc = runtmp.loc;      % locations
len = runtmp.len;      % effective length
str = runtmp.str;      % strength (kappa)
dipl_n = runtmp.did;   % diple field index

% Iterations
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% Envlope-array(x,y), Kappa-array(KX,KY), distance-array(d)
x = zeros(1,nsteps);
y = x;
d = [0:nsteps-1]*ds + min_d;


% Kappa-array
KX = zeros(1,nsteps); KY = KX;
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>nsteps )
      d2 = nsteps ;
   end;
   if ele(i)=='S'
      KX( d1:d2 ) = str(i);
      KY( d1:d2 ) = str(i);
   elseif ele(i)=='Q'
      KX( d1:d2 ) = str(i);
      KY( d1:d2 ) = -str(i);
   elseif ele(i)=='D'
      KX( d1:d2 ) = str(i)*(1-dipl_n(i));
      KY( d1:d2 ) = str(i)*dipl_n(i);
   end;
end;

X = [x0 y0 xp0 yp0];
options = optimset('LargeScale', 'on', ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'periodfunc', X,[],[],options );

xw = 1;
yw = 1;
xpw = 1;
ypw = 1;
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

