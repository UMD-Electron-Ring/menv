function [x,y,xp,yp,d,tunex,tuney] = runmenv( paramfile )
global kx ky K Ex Ey ds;
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

% Numerical parameters
max_d = runtmp.distance;
min_d = 0.0;
ds = runtmp.stepsize;          % Step-size
n = round((max_d-min_d)/ds)+1;  % steps

% Lattice
ele = runtmp.ele;      % element: 'Q'/'S'
loc = runtmp.loc;      % locations
len = runtmp.len;      % effective length
str = runtmp.str;      % strength (kappa)
dipl_n = runtmp.did;   % diple field index

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

% figure;
% plot(d,KX); hold on; plot(d,KY,'r'); hold off;
% axis([ min_d max_d min([KX,KY])*1.2 max([KX,KY])*1.2 ]);

% Leap-frog by half step
kx = KX(1); ky = KY(1);
[xpp,ypp] = calc_prim2(x0,y0);
xp(1) = xp0+xpp*ds/2;
yp(1) = yp0+ypp*ds/2;
x(1) = x0;
y(1) = y0;

% Steps
for i=1:n-1
   kx = KX(i+1); ky = KY(i+1);
   [x(i+1),y(i+1),xp(i+1),yp(i+1)] = step(x(i),y(i),xp(i),yp(i));
end;

% calculate tunesn
betax = x.^2/Ex; 
betay = y.^2/Ey;
tunex = sum(1./betax)*ds/(2*pi);
tuney = sum(1./betay)*ds/(2*pi);



% add text -- tunes
% axesHandle = findobj( gcf, 'Type', 'axes' );
% axdata = get( axesHandle(1), 'UserData' );
% if( axdata.handle(5)~=0 ) delete(axdata.handle(5)); end;
% if( axdata.handle(6)~=0 ) delete(axdata.handle(6)); end;
% yl = ylim();
% t1 = text(mean(d),0.2*(yl(2)-yl(1)),sprintf('Tune X = %.3f',tunex));
% t2 = text(mean(d),0.1*(yl(2)-yl(1)),sprintf('Tune Y = %.3f',tuney));
% axdata.handle(5)=t1; axdata.handle(6)=t2;
% set( axesHandle, 'UserData', axdata );
