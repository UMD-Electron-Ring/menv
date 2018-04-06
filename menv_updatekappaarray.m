function menv_updatekappaarray(X)
% This function updates the stored lattice model with new values according
% to an optimization routine.
%

lattmp = [];
load 'lattmp'
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1; 
loc2 = lattmp.loc2; 
OPT_ELE = lattmp.OPT_ELE;

% fudge factors (base don bench-marking with WARP bgrd elements
% Qffx = 0.955;
% Qffy = 0.935;
% Dffx = 1.9687;
% Dffy = 1;

Qffx = 1;
Qffy = 1;
Dffx = 1;
Dffy = 1;


% Update kappa
for i=1:length( X )
    if OPT_ELE(i)=='S'
        KX( loc1(i):loc2(i) ) = X(i); %0.96891*X(i);
        KY( loc1(i):loc2(i) ) = X(i); %-X(i);
    elseif OPT_ELE(i)=='Q'
        KX( loc1(i):loc2(i) ) = Qffx*X(i);
        KY( loc1(i):loc2(i) ) = -Qffy*X(i);
    elseif ele(i)=='D'
        KX( loc1(i):loc2(i) ) = Dffx*X(i)*(1-dipl_n(i));
        KY( loc1(i):loc2(i) ) = X(i)*dipl_n(i);
    end
end

lattmp.KX = KX;
lattmp.KY = KY;
save 'lattmp' lattmp

clear KX KY loc1 loc2 OPT_ELE

end
