clear all;
global Kqxx Kqyy Kqxy Kqyx K ds;

K = 1.493e-3*24/100;  % Pervence
emit = 30e-6;

% Initail conditions
x0 = 0.50414e-2;
y0 = 0.50414e-2;
xp0 = -0.022166;
yp0 = +0.022166;

% Initial 2nd order moments
xx   = (x0/2)^2;
yy   = (y0/2)^2;
xpxp = (xp0/2)^2;
ypyp = (yp0/2)^2;
xxp  = -(xx*xpxp - (emit/4)^2)^0.5;
yyp  = +(yy*ypyp - (emit/4)^2)^0.5;
xy   = 0;
xyp  = 0;
xpy  = 0;
xpyp = 0;


% Numerical parameters
rot_quad1 = 3./180*pi;
z0_quad1 = 0.08;
str_quad = 229.6;
dist_2quads = 0.16;
qlen = 0.0364;
n_quads = 16;
loc = z0_quad1:dist_2quads:(n_quads-1)*dist_2quads+z0_quad1;
len = zeros(size(loc));
len(1:end) = qlen;
str = zeros(size(loc));
str(1:2:end) = -str_quad;
str(2:2:end) = str_quad;
z_min = 0;
z_max = loc(end)+dist_2quads/2;
ds = 0.05e-3;
n = round((z_max-z_min)/ds)+1;  % steps
theta = zeros(size(loc));
theta(1) = rot_quad1;


% Kappa-array
KQXX = zeros(1,n);
KQYY = zeros(1,n);
KQXY = zeros(1,n);
KQYX = zeros(1,n);
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-z_min)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-z_min)/ds ) + 1;
   if( d2>n )
      d2 = n;
   end;
   
   KQXX( d1:d2 ) = str(i) * cos(2*theta(i));
   KQXY( d1:d2 ) = str(i) * sin(2*theta(i));
   KQYX( d1:d2 ) = KQXY( d1:d2 );
   KQYY( d1:d2 ) = -KQXX( d1:d2 );
   
end;



% Steps
XX = zeros(1,n-1);
YY = zeros(1,n-1);
XY = zeros(1,n-1);
for i=1:n-1
   Kqxx = -KQXX(i);
   Kqyy = -KQYY(i);
   Kqxy = -KQXY(i);
   Kqyx = -KQYX(i);
   [xx,xxp,xpxp,yy,yyp,ypyp,xy,xpy,xyp,xpyp] = skewstep(xx,xxp,xpxp,yy,yyp,ypyp,xy,xpy,xyp,xpyp);
   XX(i) = xx;
   YY(i) = yy;
   XY(i) = xy;
end;

x2rms = 2*sqrt(XX);
y2rms = 2*sqrt(YY);
alpha = 0.5*atan2( 2*XY, XX-YY );

figure; hold on;
plot( x2rms, 'b' );
plot( y2rms, 'r' );
figure;
plot( alpha );


