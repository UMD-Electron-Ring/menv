function [X,Y,XP,YP,D,DP] = menv_integrator(Ex,Ey,K,KX,KY,IRHO,ic,ds,allflag)
% Arguments:
% Ex is x-emittance, Ey is y-emittance
% KX,KY is array of linear focusing gradient for each step
% IRHO is array of inverse bending radius at each step
% ic is structure array holding initial conditions, in SI units
% ds is stepsize [unit meters]
% allflag = 1 returns data for every step. 
% allflag = 0 only returns data at end of last step


nsteps = length(KX);

x0 = ic.x0;
y0 = ic.y0;
xp0 = ic.xp0;
yp0 = ic.yp0;
D0 = ic.D0;
Dp0 = ic.Dp0;

kx=[]; ky=[]; % need to be initialized so they are available to all function workspaces
irho = [];

if allflag
    % -- integrate
    [X,Y,XP,YP,D,DP] = integrate();
else
    [X,Y,XP,YP,D,DP] = integrate_end();
end


% -- NESTED FUNCTIONS -- %

    function [X,Y,XP,YP,D,DP] = integrate()
        % function integrates envelope equations, saves data for every step
        
        % -- init arrays
        [X,Y,XP,YP,D,DP] = deal(zeros(1,nsteps));
        
        % Leap frog half step forward in momentum
        kx = KX(1); ky = KY(1); irho = IRHO(1);
        [xpp,ypp,dpp] = calc_prim2(x0,y0,D0);
        xp = xp0+xpp*ds/2;
        yp = yp0+ypp*ds/2;
        dp = Dp0+dpp*ds/2;
        x = x0;
        y = y0;
        d = D0;
        X(1) = x0;
        Y(1) = y0;
        D(1) = D0;
        XP(1) = xp;
        YP(1) = yp;
        DP(1) = dp;
        
        % Steps
        for i=1:nsteps-1
            kx = KX(i+1); ky = KY(i+1); irho = IRHO(i+1);
            [x,y,xp,yp,d,dp] = step(x,y,xp,yp,d,dp);
            X(i+1) = x;
            Y(i+1) = y;
            XP(i+1) = xp;
            YP(i+1) = yp;
            D(i+1) = d;
            DP(i+1) = dp;
        end
        
        % Leap frog back half step in momentum
        [xpp,ypp,dpp] = calc_prim2(x,y,d);
        xp = xp - xpp*ds/2;
        yp = yp - ypp*ds/2;
        dp = dp - dpp*ds/2;
        XP(i+1) = xp;
        YP(i+1) = yp;
        DP(i+1) = dp;
        
    end

    function [X,Y,XP,YP,D,DP] = integrate_end()
        % function integrates envelope equations, only saves/returns data for
        % last step (2x faster than saving all steps)
        
        
        % Leap frog half step
        kx = KX(1); ky = KY(1); irho = IRHO(1);
        [xpp,ypp,dpp] = calc_prim2(x0,y0,D0);
        xp = xp0+xpp*ds/2;
        yp = yp0+ypp*ds/2;
        dp = Dp0+dpp*ds/2;
        x = x0;
        y = y0;
        d = D0;
        
        % Steps
        for i=1:nsteps-1
            kx = KX(i+1); ky = KY(i+1); irho = IRHO(1+1);
            [x,y,xp,yp,d,dp] = step(x,y,xp,yp,d,dp);
        end
        
        % xp and yp back half step
        [xpp,ypp,dpp] = calc_prim2(x,y,d);
        xp = xp - xpp*ds/2;
        yp = yp - ypp*ds/2;
        dp = dp - dpp*ds/2;
        
        % return last step
        X = x; Y = y; XP = xp; YP = yp; D = d; DP = dp;
        
    end

    function [x2,y2,xp2,yp2,d2,dp2] = step(x1,y1,xp1,yp1,d1,dp1)
        % function for single step
        x2 = x1+xp1*ds;
        y2 = y1+yp1*ds;
        d2 = d1+dp1*ds;
        [xpp,ypp,dpp] = calc_prim2(x2,y2,d2);
        xp2 = xp1+xpp*ds;
        yp2 = yp1+ypp*ds;
        dp2 = dp1+dpp*ds;
    end

    function [xpp1,ypp1,dpp1] = calc_prim2(x1,y1,d1)
        % function to calculate velocity impulse
        % Reiser eq. 4.191, 4.192
        xpp1 = -( kx*x1-2*K/(x1+y1)-Ex^2/(x1^3) );
        ypp1 = -( ky*y1-2*K/(x1+y1)-Ey^2/(y1^3) );
        dpp1 = irho - kx*d1;
    end
  
end
