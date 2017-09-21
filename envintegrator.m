function [X,Y,XP,YP] = envintegrator(KX,KY,ic,ds,nsteps,allflag)

global clm
K = clm.usrdata.perveance; % Pervence
Ex = clm.usrdata.emitance; % Emmitance x
Ey = Ex;

x0 = ic(1);
y0 = ic(2);
xp0 = ic(3);
yp0 = ic(4);

kx=[]; ky=[]; % need to be initialized so they are available to all function workspaces

if allflag
    % -- init arrays
    [X,Y,XP,YP] = deal(zeros(1,nsteps));
    X(1) = x0;
    Y(1) = y0;
    % -- integrate
    [X,Y,XP,YP] = integrate();
else
    [X,Y,XP,YP] = integrate_end();
end


% -- NESTED FUNCTIONS -- %

    function [X,Y,XP,YP] = integrate()
        % function integrates envelope equations, saves data for every step
        
        % Leap frog half step
        kx = KX(1); ky = KY(1);
        [xpp,ypp] = calc_prim2(x0,y0);
        xp = xp0+xpp*ds/2;
        yp = yp0+ypp*ds/2;
        x = x0;
        y = y0;
        XP(1) = xp;
        YP(1) = yp;
        
        % Steps
        for i=1:nsteps-1
            kx = KX(i+1); ky = KY(i+1);
            [x,y,xp,yp] = step(x,y,xp,yp);
            X(i+1) = x;
            Y(i+1) = y;
            XP(i+1) = xp;
            YP(i+1) = yp;
        end
        
    end

    function [X,Y,XP,YP] = integrate_end()
        % function integrates envelope equations, only saves/returns data for
        % last step (2x faster than saving all steps)
        
        
        % Leap frog half step
        kx = KX(1); ky = KY(1);
        [xpp,ypp] = calc_prim2(x0,y0);
        xp = xp0+xpp*ds/2;
        yp = yp0+ypp*ds/2;
        x = x0;
        y = y0;
        
        % Steps
        for i=1:nsteps-1
            kx = KX(i+1); ky = KY(i+1);
            [x,y,xp,yp] = step(x,y,xp,yp);
        end
        
        % xp and yp back half step
        [xpp,ypp] = calc_prim2(x,y);
        xp = xp - xpp*ds/2;
        yp = yp - ypp*ds/2;
        
        % return last step
        X = x; Y = y; XP = xp; YP = xp;
        
    end

    function [x,y,xp,yp] = step(x0,y0,xp0,yp0)
        % function for single step
        x = x0+xp0*ds;
        y = y0+yp0*ds;
        [xpp,ypp] = calc_prim2(x,y);
        xp = xp0+xpp*ds;
        yp = yp0+ypp*ds;
    end

    function [xpp,ypp] = calc_prim2(x,y)
        % function to calculate velocity impulse
        xpp = -( kx*x-2*K/(x+y)-Ex^2/(x^3) );
        ypp = -( ky*y-2*K/(x+y)-Ey^2/(y^3) );
    end
  
end
