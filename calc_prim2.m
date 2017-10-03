function [xpp,ypp,Dpp] = calc_prim2(x,y,D)
global kx ky K Ex Ey irho;
xpp = -( kx*x-2*K/(x+y)-Ex^2/(x^3) );
ypp = -( ky*y-2*K/(x+y)-Ey^2/(y^3) );
Dpp = irho - kx*D;
