function [xpp,ypp] = calc_prim2(x,y)
global kx ky K Ex Ey;
xpp = -( kx*x-2*K/(x+y)-Ex^2/(x^3) );
ypp = -( ky*y-2*K/(x+y)-Ey^2/(y^3) );
