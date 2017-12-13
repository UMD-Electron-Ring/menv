function X = match2target( paramfile )
global kx ky irho K Ex Ey ds;
global loc1 loc2 KX KY IRHO nsteps;
global x0 y0 xp0 yp0 D0 Dp0;
global x1 y1 xp1 yp1 D1 Dp1;
global xw yw xpw ypw Dw Dpw;
global OPT_ELE;
load 'runtmp';

% Global beam parateters
K = runtmp.perveance; % Pervence
Ex = runtmp.emitance; % Emmitance x
Ey = Ex;               % Emmitance y
% Initial conditions
x0 = runtmp.x0;
y0 = runtmp.y0;
xp0 = runtmp.xp0;
yp0 = runtmp.yp0;
D0 = runtmp.D0;
Dp0 = runtmp.Dp0;
% target conditions
x1 = runtmp.x1;
y1 = runtmp.y1;
xp1 = runtmp.xp1;
yp1 = runtmp.yp1;
D1 = runtmp.D1;
Dp1 = runtmp.Dp1;
% weights
xw = runtmp.xw;
yw = runtmp.yw;
xpw = runtmp.xpw;
ypw = runtmp.ypw;
Dw = runtmp.Dw;
Dpw = runtmp.Dpw;
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
invrho = runtmp.irho;

% Optimize which
opt = runtmp.opt;
% Iterations
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% Envlope-array(x,y), Kappa-array(KX,KY), distance-array(d)
[x,y,D] = deal(zeros(1,nsteps));
d = [0:nsteps-1]*ds + min_d;


% Kappa-array
[KX,KY,IRHO] = deal(zeros(1,nsteps)); 
loc1 = []; loc2 = []; X0 = []; OPT_ELE = [];
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>nsteps )
      d2 = nsteps ;
   end;
   if( opt(i)==0 )
      if ele(i)=='S'
         KX( d1:d2 ) = str(i);
         KY( d1:d2 ) = str(i);
         IRHO( d1:d2 ) = invrho(i);
      elseif ele(i)=='Q'
         KX( d1:d2 ) = 0.955*str(i);
         KY( d1:d2 ) = -0.935*str(i);
         IRHO( d1:d2 ) = invrho(i);
      elseif ele(i)=='D'
         KX( d1:d2 ) = 1.9687*str(i)*(1-dipl_n(i));
         KY( d1:d2 ) = str(i)*dipl_n(i);
         IRHO( d1:d2 ) = invrho(i);
      end;
   else
      loc1 = [ loc1, d1 ];
      loc2 = [ loc2, d2 ];
      X0 = [ X0, str(i) ];
      OPT_ELE = [ OPT_ELE, ele(i) ];
   end;
end;

X = X0;
if( length(X)<=4 ) scale = 'on';
else scale = 'off'; end;
options = optimset('LargeScale', scale, ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'stepfunc', X,[],[],options );


f =  stepfunc( X ); % run just to get error contributions;
runtmp.f = f;
save 'runtmp' runtmp;

