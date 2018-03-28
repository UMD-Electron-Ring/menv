function [x,y,xp,yp,D,Dp,d,tunex,tuney,Cx,Cy] = runmenv()

% -- make lattmp variable w/ lattice arrays
X0 = menv_makekappaarray(); clear X0;

% -- load lattice model parameters
load 'lattmp'
IRHO = lattmp.IRHO;
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1;
loc2 = lattmp.loc2;
OPT_ELE = lattmp.OPT_ELE;
nsteps = lattmp.nsteps;
d = lattmp.d;

% -- load beam parameters from rntump
load 'runtmp'
K = runtmp.perveance; % Perveance
Ex = runtmp.emitance; % Emmitance x
Ey = Ex;
ds = runtmp.stepsize;

% -- integrate envelope
ic = runtmp.ic;
[x,y,xp,yp,D,Dp] = menv_integrator(Ex,Ey,K,KX,KY,IRHO,ic,ds,1);


% calculate tunes
betax = x.^2/Ex; 
betay = y.^2/Ey;
tunex = sum(1./betax)*ds/(2*pi);
tuney = sum(1./betay)*ds/(2*pi);

% Calculate chromaticity
Cx = -1/(4*pi)*sum(KX.*betax)*ds;
Cy = -1/(4*pi)*sum(KY.*betay)*ds;



% add text to figure -- display tunes
% axesHandle = findobj( gcf, 'Type', 'axes' );
% axdata = get( axesHandle(1), 'UserData' );
% if( axdata.handle(5)~=0 ) delete(axdata.handle(5)); end;
% if( axdata.handle(6)~=0 ) delete(axdata.handle(6)); end;
% yl = ylim();
% t1 = text(mean(d),0.2*(yl(2)-yl(1)),sprintf('Tune X = %.3f',tunex));
% t2 = text(mean(d),0.1*(yl(2)-yl(1)),sprintf('Tune Y = %.3f',tuney));
% axdata.handle(5)=t1; axdata.handle(6)=t2;
% set( axesHandle, 'UserData', axdata );
