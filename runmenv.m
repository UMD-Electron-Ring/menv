function [x,y,d] = runmenv( paramfile )
global kx ky K Ex Ey ds;
load( paramfile );

% Global beam parateters
K = usrdata.perveance; % Pervence
Ex = usrdata.emitance; % Emmitance x
Ey = Ex;               % Emmitance y

% Initail conditions
x0 = usrdata.x0;
y0 = usrdata.y0;
xp0 = usrdata.xp0;
yp0 = usrdata.yp0;

% Numerical parameters
max_d = usrdata.distance;
min_d = 0.0;
ds = usrdata.stepsize;          % Step-size
n = round((max_d-min_d)/ds)+1;  % steps

% Lattice
ele = usrdata.ele;      % element: 'Q'/'S'
loc = usrdata.loc;      % locations
len = usrdata.len;      % effective length
str = usrdata.str;      % strength (kappa)
dipl_n = usrdata.did;   % diple field index

% Envlope-array(x,y), Kappa-array(KX,KY), distance-array(d)
x = zeros(1,n);
y = x;
d = [0:n-1]*ds + min_d;

% Kappa-array
KX = zeros(1,n); KY = KX;
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>n )
      d2 = n;
   end;
   if ele(i)=='S'
      KX( d1:d2 ) = str(i); %str(i)*0.96891;
      KY( d1:d2 ) = str(i); %-str(i);
   elseif ele(i)=='Q'
      KX( d1:d2 ) = str(i);
      KY( d1:d2 ) = -str(i);
   elseif ele(i)=='D'
      KX( d1:d2 ) = str(i)*(1-dipl_n(i));
      KY( d1:d2 ) = str(i)*dipl_n(i);
   end;
end;

%figure;
%plot(d,KX); hold on; plot(d,KY,'r'); hold off;
%axis([ min_d max_d min([KX,KY])*1.2 max([KX,KY])*1.2 ]);

% Leap-frog by half step
kx = KX(1); ky = KY(1);
[xpp,ypp] = calc_prim2(x0,y0);
xp = xp0+xpp*ds/2;
yp = yp0+ypp*ds/2;
x(1) = x0;
y(1) = y0;

% Steps
for i=1:n-1
   kx = KX(i+1); ky = KY(i+1);
   [x(i+1),y(i+1),xp,yp] = step(x(i),y(i),xp,yp);
end;

% plot envelope
%figure( fig );
%plot(d,x); hold on; plot(d,y,'r'); hold off;
%axis([ min_d max_d 0.0 max([x,y])*1.2 ]);
