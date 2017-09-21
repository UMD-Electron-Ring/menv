function f =  stepfunc2( X )
global KX KY
global loc1 loc2 
global OPT_ELE;
global clm

x0 = clm.usrdata.x0;
y0 = clm.usrdata.y0;
xp0 = clm.usrdata.xp0;
yp0 = clm.usrdata.yp0;
x1 = clm.usrdata.x1;
y1 = clm.usrdata.y1;
xp1 = clm.usrdata.xp1;
yp1 = clm.usrdata.yp1;
xw = clm.usrdata.xw;
yw = clm.usrdata.yw;
xpw = clm.usrdata.xpw;
ypw = clm.usrdata.ypw;
emittance = clm.usrdata.emitance;

nuxt = clm.usrdata.nuxt;
nuyt = clm.usrdata.nuyt;
nuxw = clm.usrdata.nuxw;
nuyw = clm.usrdata.nuyw;
betaw = clm.usrdata.betaw;

ds = clm.usrdata.stepsize;
distance = clm.usrdata.distance;
stepsize = clm.usrdata.stepsize;
nsteps = round((distance)/ds)+1;  % steps

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

% -- integrate
[x,y,xp,yp] = envintegrator(KX,KY,[x0,y0,xp0,yp0],ds,nsteps,1);

% -- calculate tune values
Lchan = 32;
d = 0:ds:distance;

betax = (x).^2/(emittance);
betay = (y).^2/(emittance);

% -- full lattice tune
psix = stepsize*sum(1./betax);
psiy = stepsize*sum(1./betay);
nux = psix/(2*pi);
nuy = psiy/(2*pi);

% condition for integer ring tune
nuxcontr = nux-nuxt;
nuycontr = nuy-nuyt;

% condition for envelope target
xcontr = x(end)-x1;
ycontr = y(end)-y1;
xpcontr = xp(end)-xp1;
ypcontr = yp(end)-yp1; 

% condition for symmetric beam
betacontr = max(abs(betax-betay));

f = [xcontr,ycontr,xpcontr,ypcontr,nuxcontr,nuycontr,betacontr];
f = f.*[xw yw xpw ypw nuxw nuyw betaw];
