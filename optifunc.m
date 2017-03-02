function f =  optifunc( X )
global Ex Ey kx ky ds;
global loc1 loc2 KX KY nsteps;
global x0 y0 xp0 yp0;
global x1 y1 xp1 yp1;
global xw yw xpw ypw;
global yref refw tunew;
global OPT_ELE;


% Evaluate kappa
for i=1:length( X )
    if OPT_ELE(i)=='S'
        KX( loc1(i):loc2(i) ) = X(i); %0.96891*X(i);
        KY( loc1(i):loc2(i) ) = X(i); %-X(i);
    elseif OPT_ELE(i)=='Q'
        KX( loc1(i):loc2(i) ) = X(i);
        KY( loc1(i):loc2(i) ) = -X(i);
    elseif ele(i)=='D'
        KX( loc1(i):loc2(i) ) = X(i)*(1-dipl_n(i));
        KY( loc1(i):loc2(i) ) = X(i)*dipl_n(i);
    end;
end;

% allocate history arrays
[x,y,xp,yp] = deal(zeros(1,nsteps));
d = [0:nsteps-1]*ds;


% Leap frog half step
kx = KX(1); ky = KY(1);
[xpp,ypp] = calc_prim2(x0,y0);
xp(1) = xp0+xpp*ds/2;
yp(1) = yp0+ypp*ds/2;
x(1) = x0;
y(1) = y0;

% Steps
for i=1:nsteps-1
   kx = KX(i+1); ky = KY(i+1);
   [x(i+1),y(i+1),xp(i+1),yp(i+1)] = step(x(i),y(i),xp(i),yp(i));
end

% xp and yp back half step
[xpp,ypp] = calc_prim2(x(end),y(end));
xp(end) = xp(end-1) - xpp*ds/2;
yp(end) = yp(end-1) - ypp*ds/2;


% calculate tunes
betax = x.^2/Ex; 
betay = y.^2/Ey;
tunex = sum(1./betax)*ds/(2*pi);
tuney = sum(1./betay)*ds/(2*pi);

ind = find(abs(d-.32)==min(abs(d-.32)));
tunex_dr = sum(1./betax(1:ind))*ds/(2*pi);
tuney_dr = sum(1./betay(1:ind))*ds/(2*pi);
tunex_ring = 8*tunex-2*tunex_dr;
tuney_ring = 8*tuney-2*tuney_dr;
tunex_res = mod(tunex_ring,1);
tuney_res = mod(tuney_ring,1);


% minimization function
f = [x(end)-x1,y(end)-y1,xp(end)-xp1,yp(end)-yp1,mean(abs(x-yref)),mean(abs(y-yref)),tunex-.1522,tuney-.1525];
f = f.*[xw yw xpw ypw refw refw tunew tunew];



% plot
x=x*1.e2; y=y*1.e2; d=d*1.e2;  % m-cm
axesHandle = findobj( gcf, 'Type', 'axes' );
axdata = get( axesHandle(1), 'UserData' );
if( axdata.handle(3)~=0 ) delete(axdata.handle(3)); end;
if( axdata.handle(4)~=0 ) delete(axdata.handle(4)); end;
hold on; h3 = plot(d,x,':b'); h4 = plot(d,y,':r'); hold off;
drawnow();
axdata.handle(3)=h3; axdata.handle(4)=h4; 

% add text -- tunes
if( axdata.handle(5)~=0 ) delete(axdata.handle(5)); end;
if( axdata.handle(6)~=0 ) delete(axdata.handle(6)); end;
yl = ylim();
t1 = text(mean(d),0.2*(yl(2)-yl(1)),sprintf('Tune X = %.3f',tunex));
t2 = text(mean(d),0.1*(yl(2)-yl(1)),sprintf('Tune Y = %.3f',tuney));
axdata.handle(5)=t1; axdata.handle(6)=t2;
set( axesHandle, 'UserData', axdata );

end
