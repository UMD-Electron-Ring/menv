function menv_updatekappaarray(X)
% This function updates the stored lattice model with new values according
% to an optimization routine.
%

% -- this is a stop-gap measure until I figure out what better way to
% manage passing structure variables without hard-disk writes
global lattmp 
%lattmp = [];
%load 'lattmp'
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1; 
loc2 = lattmp.loc2; 
OPT_ELE = lattmp.OPT_ELE;

% fudge factors (based on bench-marking with WARP bgrd elements)
% Qffx = 0.955;
% Qffy = 0.935;
% Dffx = 1.9687;
% Dffy = 1;

Qffx = 1;
Qffy = 1;
Dffx = .5*.146/.556; % -- horz focusing is weaker than vert. focusing. From RK dipole note + 2010 online table
Dffy = 1;
PDffx = 0.00/0.0410; % -- horz focusing is weaker than vert. focusing. From RK dipole note + 2010 online table
PDffy = 1;


% Update kappa
for i=1:length( X )
    if OPT_ELE(i)=='S' % -- solenoid
        KX( loc1(i):loc2(i) ) = X(i); 
        KY( loc1(i):loc2(i) ) = X(i);
    elseif OPT_ELE(i)=='Q'  % -- ring quad
        KX( loc1(i):loc2(i) ) = Qffx*X(i);
        KY( loc1(i):loc2(i) ) = -Qffy*X(i);
    elseif ele(i)=='D' % -- ring dipole
        KX( loc1(i):loc2(i) ) = Dffx*X(i);
        KY( loc1(i):loc2(i) ) = Dffy*X(i);
    elseif OPT_ELE(i)=='Y' % -- enlarged Y-section quad
        KX( loc1(i):loc2(i) ) = X(i);
        KY( loc1(i):loc2(i) ) = -X(i);
    elseif ele(i)=='P' % -- enlarged Y-section dipole
        KX( loc1(i):loc2(i) ) = PDffx*X(i);
        KY( loc1(i):loc2(i) ) = PDffy*X(i);
    end
end

lattmp.KX = KX;
lattmp.KY = KY;
%save('lattmp.mat','lattmp')

clear KX KY loc1 loc2 OPT_ELE

end
