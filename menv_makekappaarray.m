function [X0] = menv_makekappaarray()
% this function makes a lattice model with the resolution determined by
% stepsize, using the information in the lattice definition in clm.usrdata
% It returns a list of focusing strengths (X0) to be used for optimization

% Load numerical parameters
load( 'runtmp' );
max_d = runtmp.distance;
min_d = 0.0;
ds = runtmp.stepsize;          % Step-size
nsteps = round((max_d-min_d)/ds)+1;  % steps
% Load Lattice parameters
ele = runtmp.ele;      % element: 'Q'/'S'
loc = runtmp.loc;      % locations
len = runtmp.len;      % effective length
str = runtmp.str;      % strength (kappa)
dipl_n = runtmp.did;   % diple field index
try % -- irho is not defined for old .spt files, if undefined use 0-array 
    irho = runtmp.irho;    % inverse bend radius
catch
    irho = runtmp.str*0;
end
opt = runtmp.opt;      % flag to use element in optimization


% distance-array(d)
d = [0:nsteps-1]*ds + min_d;

% Kappa-array
[KX,KY,IRHO] = deal(zeros(1,nsteps)); 
loc1 = []; loc2 = []; X0 = []; OPT_ELE = [];
for i=1:length(loc)
   d1 = round( (loc(i)-len(i)/2-min_d)/ds ) + 1;
   d2 = round( (loc(i)+len(i)/2-min_d)/ds ) + 1;
   if( d2>nsteps )
      d2 = nsteps ;
   end
   if( opt(i)==0 )
      if ele(i)=='S'
         KX( d1:d2 ) = str(i);
         KY( d1:d2 ) = str(i);
         IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='Q'
         KX( d1:d2 ) = 0.955*str(i);
         KY( d1:d2 ) = -0.935*str(i);
         IRHO( d1:d2 ) = irho(i);
      elseif ele(i)=='D'
         KX( d1:d2 ) = 1.9687*str(i)*(1-dipl_n(i));
         KY( d1:d2 ) = str(i)*dipl_n(i);
         IRHO( d1:d2 ) = irho(i);
      end
   else
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
save 'lattmp' lattmp