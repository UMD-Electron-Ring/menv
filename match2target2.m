function X = match2target2()
global loc1 loc2 KX KY;
global OPT_ELE;
global clm

% Global beam parateters
K = clm.usrdata.perveance; % Pervence
Ex = clm.usrdata.emitance; % Emmitance x
Ey = Ex;               % Emmitance y
% Initial conditions
x0 = clm.usrdata.x0;
y0 = clm.usrdata.y0;
xp0 = clm.usrdata.xp0;
yp0 = clm.usrdata.yp0;
% target conditions
x1 = clm.usrdata.x1;
y1 = clm.usrdata.y1;
xp1 = clm.usrdata.xp1;
yp1 = clm.usrdata.yp1;
% weights
xw = clm.usrdata.xw;
yw = clm.usrdata.yw;
xpw = clm.usrdata.xpw;
ypw = clm.usrdata.ypw;
% Numerical parameters
max_d = clm.usrdata.distance;
min_d = 0.0;
ds = clm.usrdata.stepsize;          % Step-size
nsteps = round((max_d-min_d)/ds)+1;  % steps
% Lattice
ele = clm.usrdata.ele;      % element: 'Q'/'S'
loc = clm.usrdata.loc;      % locations
len = clm.usrdata.len;      % effective length
str = clm.usrdata.str;      % strength (kappa)
dipl_n = clm.usrdata.did;   % diple field index

% Optimize which
opt = clm.usrdata.opt;
% Iterations
maxIter = round(clm.usrdata.maxIter);
tolFun = clm.usrdata.tolFun;

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

X = X0; % X is values of optimization variables
if( length(X)<=4 ) scale = 'on';
else scale = 'off'; end;
options = optimset('LargeScale', scale, ...
   'Display', 'iter', ...
   'MaxIter', maxIter, ...
   'TolFun', tolFun );
X = lsqnonlin( 'stepfunc2', X,[],[],options );

