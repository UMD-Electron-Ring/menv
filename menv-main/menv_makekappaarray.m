function [X0] = menv_makekappaarray()
% this function makes a lattice model with the resolution determined by
% stepsize, using the information in the lattice definition in clm.usrdata
% It returns a list of focusing strengths (X0) to be used for optimization

% Load numerical parameters
% -- this is a stop-gap measure until I figure out what better way to
% manage passing structure variables without hard-disk writes
% load( 'runtmp' );
global lattmp runtmp

max_d = runtmp.distance;
try min_d = runtmp.d0;
catch min_d = 0.0;
end
ds = runtmp.stepsize;          % Step-size
nsteps = round((max_d-min_d)/ds)+1;  % steps
% Load Lattice parameters
ele = runtmp.ele;      % element: 'Q'/'S'
loc = runtmp.loc;      % locations
len = runtmp.len;      % effective length
str = runtmp.str;      % strength (kappa)
try % -- irho is not defined for old .spt files, if undefined use 0-array 
    irho = runtmp.irho;    % inverse bend radius
catch
    irho = runtmp.str*0;
end
opt = runtmp.opt;      % flag to use element in optimization

% distance-array(d)
d = [0:nsteps-1]*ds + min_d;

Qffx = 1;
Qffy = 1;
Dffx = 0; % -- horz defocusing canceled by geom. focusing
%Dffx = .5*.146/.556; % -- horz defocusing is weaker than vert. focusing. From RK dipole note + 2010 online table
Dffy = 1;
PDffx = 0.00/0.0410; % -- horz defocusing is weaker than vert. focusing. From RK dipole note + 2010 online table
PDffy = 1;


% Kappa-array 
% includes fudge-factors, which should be moved to a different place
% (current2kappa probably)
[KX,KY,IRHO] = deal(zeros(1,nsteps)); 
loc1 = []; loc2 = []; X0 = []; OPT_ELE = [];
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>nsteps )
      d2 = nsteps ;
   end
      if ele(i)=='S' % -- solenoid
         KX( d1:d2 ) = str(i);
         KY( d1:d2 ) = str(i);
         IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='Q' % -- ring quad
         KX( d1:d2 ) = Qffx*str(i);
         KY( d1:d2 ) = -Qffy*str(i);
         IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='D' % -- ring dipole
         KX( d1:d2 ) = Dffx*str(i);
         KY( d1:d2 ) = Dffy*str(i);
         IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='Y' % -- enlarged Y-section quad
          KX( d1:d2 ) = str(i);
          KY( d1:d2 ) = -str(i);
          IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='P' % -- enlarged Y-section dipole
          KX( d1:d2 ) = PDffx*str(i);
          KY( d1:d2 ) = PDffy*str(i);
          IRHO( d1:d2 ) = irho(i);
      end
   if( opt(i)==1 )
      loc1 = [ loc1, d1 ];
      loc2 = [ loc2, d2 ];
      X0 = [ X0, str(i) ];
      OPT_ELE = [ OPT_ELE, ele(i) ];
   end
end


% -- save data to file
lattmp = struct();
lattmp.OPT_ELE = OPT_ELE;
lattmp.loc1 = loc1;
lattmp.loc2 = loc2;
lattmp.IRHO = IRHO;
lattmp.nsteps = nsteps;
lattmp.KX = KX;
lattmp.KY = KY;
lattmp.d = d;
%save('lattmp.mat','lattmp')