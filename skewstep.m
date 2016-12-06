function [xx1,xxp1,xpxp1,yy1,yyp1,ypyp1,xy1,xpy1,xyp1,xpyp1] = skewstep(xx,xxp,xpxp,yy,yyp,ypyp,xy,xpy,xyp,xpyp)
global Kqxx Kqyy Kqxy Kqyx K ds;

a = 0.5 * atan2( 2*xy, xx-yy );
xbxb = xx*cos(a)^2 + yy*sin(a)^2 + 2*xy*cos(a)*sin(a);
ybyb = yy*cos(a)^2 + xx*sin(a)^2 - 2*xy*cos(a)*sin(a);
Ksxb = K / 2 / (xbxb + (xbxb*ybyb)^0.5);
Ksyb = K / 2 / (ybyb + (xbxb*ybyb)^0.5);

Ksxx = Ksxb*cos(a)^2 + Ksyb*sin(a)^2;
Ksxy = (Ksxb-Ksyb) * sin(a)*cos(a);
Ksyx = Ksxy;
Ksyy = Ksyb*cos(a)^2 + Ksxb*sin(a)^2;

Kxx = Kqxx + Ksxx;
Kyy = Kqyy + Ksyy;
Kxy = Kqxy + Ksxy;
Kyx = Kxy;

xx1   = xx   + ds * (2*xxp);
xxp1  = xxp  + ds * (xpxp + Kxx*xx + Kxy*xy);
xpxp1 = xpxp + ds * (2*Kxx*xxp + 2*Kxy*xpy);

yy1   = yy   + ds * (2*yyp);
yyp1  = yyp  + ds * (ypyp + Kyy*yy + Kyx*xy);
ypyp1 = ypyp + ds * (2*Kyy*yyp + 2*Kyx*xyp);

xy1   = xy   + ds * (xpy + xyp);
xpy1  = xpy  + ds * (xpyp + Kxx*xy + Kxy*yy);
xyp1  = xyp  + ds * (xpyp + Kyy*xy + Kyx*xx);
xpyp1 = xpyp + ds * (Kxx*xyp + Kxy*yyp + Kyy*xpy + Kyx*xxp);
