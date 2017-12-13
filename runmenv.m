function [x,y,xp,yp,D,Dp,d,tunex,tuney,Cx,Cy] = runmenv( paramfile )
global kx ky K Ex Ey ds irho;
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
D0 = runtmp.D0;
Dp0 = runtmp.Dp0;

% Numerical parameters
max_d = runtmp.distance;
min_d = runtmp.s0;
ds = runtmp.stepsize;          % Step-size
n = round((max_d-min_d)/ds)+1;  % steps

% Lattice
ele = runtmp.ele;      % element: 'Q'/'S'
loc = runtmp.loc;      % locations
len = runtmp.len;      % effective length
str = runtmp.str;      % strength (kappa)
dipl_n = runtmp.did;   % diple field index
invrho = runtmp.irho; % inverse rho

% Envlope-array(x,y), Kappa-array(KX,KY), distance-array(d)
x = zeros(1,n);
y = x;
d = [0:n-1]*ds + min_d;

% Kappa-array
[KX,KY,IRHO] = deal(zeros(1,n));
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>n )
      d2 = n;
   end;
   % -- break if d2 < 0 (before run start)
   if d2<0
       continue
   end
   % -- if starting mid-element, applying partial element as needed
   if d1<0
       d1=1
   end
   if ele(i)=='S'
       KX( d1:d2 ) = str(i); %str(i)*0.96891;
       KY( d1:d2 ) = str(i); %-str(i);
       IRHO( d1:d2 ) = invrho(i);
   elseif ele(i)=='Q'
       % -- fudge factors for agreement w/ WARP env. model
       KX( d1:d2 ) = .955*str(i);
       KY( d1:d2 ) = -.935*str(i);
       IRHO( d1:d2 ) = invrho(i);
   elseif ele(i)=='D'
       KX( d1:d2 ) = 1.9687*str(i)*(1-dipl_n(i)); % 1.9687 is fudge factor used to get agreement w/ WARP env model
       KY( d1:d2 ) = str(i)*dipl_n(i);
       IRHO( d1:d2 ) = invrho(i);
   end;
end;

% figure;
% plot(d,KX); hold on; plot(d,KY,'r'); hold off;
% axis([ min_d max_d min([KX,KY])*1.2 max([KX,KY])*1.2 ]);

% Leap-frog by half step
kx = KX(1); ky = KY(1); irho = IRHO(1);
[xpp,ypp,Dpp] = calc_prim2(x0,y0,D0);
xp(1) = xp0+xpp*ds/2;
yp(1) = yp0+ypp*ds/2;
Dp(1) = Dp0+Dpp*ds/2;
x(1) = x0;
y(1) = y0;
D(1) = D0;

% Steps
for i=1:n-1
   kx = KX(i+1); ky = KY(i+1); irho = IRHO(i+1);
   [x(i+1),y(i+1),xp(i+1),yp(i+1),D(i+1),Dp(i+1)] = step(x(i),y(i),xp(i),yp(i),D(i),Dp(i));
end;

% calculate tunesn
betax = x.^2/Ex; 
betay = y.^2/Ey;
tunex = sum(1./betax)*ds/(2*pi);
tuney = sum(1./betay)*ds/(2*pi);

% Calculate chromaticity

Cx = -1/(4*pi)*sum(KX.*betax)*ds;
Cy = -1/(4*pi)*sum(KY.*betay)*ds;



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
