function X = match2target( paramfile )
global kx ky K Ex Ey ds;
global loc1 loc2 KX KY nsteps;
global x0 y0 xp0 yp0;
global x1 y1 xp1 yp1;
global xw yw xpw ypw;
global OPT_ELE;
load( paramfile );

% Global beam parateters
K = usrdata.perveance; % Pervence
Ex = usrdata.emitance; % Emmitance x
Ey = Ex;               % Emmitance y
% Initial conditions
x0 = usrdata.x0;
y0 = usrdata.y0;
xp0 = usrdata.xp0;
yp0 = usrdata.yp0;
% target conditions
x1 = usrdata.x1;
y1 = usrdata.y1;
xp1 = usrdata.xp1;
yp1 = usrdata.yp1;
% weights
xw = usrdata.xw;
yw = usrdata.yw;
xpw = usrdata.xpw;
ypw = usrdata.ypw;
% Numerical parameters
max_d = usrdata.distance;
min_d = 0.0;
ds = usrdata.stepsize;          % Step-size
nsteps = round((max_d-min_d)/ds)+1;  % steps
% Lattice
ele = usrdata.ele;      % element: 'Q'/'S'
loc = usrdata.loc;      % locations
len = usrdata.len;      % effective length
str = usrdata.str;      % strength (kappa)
dipl_n = usrdata.did;   % diple field index

% Optimize which
opt = usrdata.opt;
% Iterations
maxIter = round(usrdata.maxIter);
tolFun = usrdata.tolFun;

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

X = X0;
if( length(X)<=4 ) scale = 'on';
else scale = 'off'; end;
options = optimset('LargeScale', scale, ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'stepfunc', X,[],[],options );

