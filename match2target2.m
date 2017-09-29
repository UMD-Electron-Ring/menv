function X = match2target2( paramfile )

global loc1 loc2 KX KY;
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
% target conditions
x1 = runtmp.x1;
y1 = runtmp.y1;
xp1 = runtmp.xp1;
yp1 = runtmp.yp1;
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

% Optimize which
opt = runtmp.opt;
% Iterations
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% Envlope-array(x,y), Kappa-array(KX,KY), distance-array(d)
x = zeros(1,nsteps);
y = x;
d = [0:nsteps-1]*ds + min_d;


% Kappa-array
KX = zeros(1,nsteps); KY = KX;
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
      elseif ele(i)=='Q'
         KX( d1:d2 ) = str(i);
         KY( d1:d2 ) = -str(i);
      elseif ele(i)=='D'
         KX( d1:d2 ) = str(i)*(1-dipl_n(i));
         KY( d1:d2 ) = str(i)*dipl_n(i);
      end;
   else
      loc1 = [ loc1, d1 ];
      loc2 = [ loc2, d2 ];
      X0 = [ X0, str(i) ];
      OPT_ELE = [ OPT_ELE, ele(i) ];
   end;
end;

runtmp.KX = KX;
runtmp.KY = KY;
runtmp.OPT_ELE = OPT_ELE;
runtmp.loc1=loc1;
runtmp.loc2=loc2;
save 'runtmp' runtmp

X = X0; % X is values of optimization variables
if( length(X)<=4 ) scale = 'on';
else scale = 'off'; end;
options = optimset('LargeScale', scale, ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'stepfunc2', X,[],[],options );

