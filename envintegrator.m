function [X,Y,XP,YP] = envintegrator(KX,KY,ic,ds,nsteps,allflag)

runtmp = []; % needed for fixed workspace
load 'runtmp'
K = runtmp.perveance; % Pervence
Ex = runtmp.emitance; % Emmitance x
Ey = Ex;

x0 = ic(1);
y0 = ic(2);
xp0 = ic(3);
yp0 = ic(4);

kx=[]; ky=[]; % need to be initialized so they are available to all function workspaces

if allflag
    % -- integrate
    [X,Y,XP,YP] = integrate();
else
    [X,Y,XP,YP] = integrate_end();
end


% -- NESTED FUNCTIONS -- %

    function [X,Y,XP,YP] = integrate()
        % function integrates envelope equations, saves data for every step
        
        % -- init arrays
        [X,Y,XP,YP] = deal(zeros(1,nsteps));
        
        % Leap frog half step
        kx = KX(1); ky = KY(1);
        [xpp,ypp] = calc_prim2(x0,y0);
        xp = xp0+xpp*ds/2;
        yp = yp0+ypp*ds/2;
        x = x0;
        y = y0;
        X(1) = x0;
        Y(1) = y0;
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

    function [x2,y2,xp2,yp2] = step(x1,y1,xp1,yp1)
        % function for single step
        x2 = x1+xp1*ds;
        y2 = y1+yp1*ds;
        [xpp,ypp] = calc_prim2(x2,y2);
        xp2 = xp1+xpp*ds;
        yp2 = yp1+ypp*ds;
    end

    function [xpp1,ypp1] = calc_prim2(x1,y1)
        % function to calculate velocity impulse
        xpp1 = -( kx*x1-2*K/(x1+y1)-Ex^2/(x1^3) );
        ypp1 = -( ky*y1-2*K/(x1+y1)-Ey^2/(y1^3) );
    end
  
end
