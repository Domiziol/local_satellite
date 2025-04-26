clear;

hm = 6378137;
hs = 20000000;
vs = 299792458;

r = hm+hs;
theta = [52.885907 50.312052 47.796902 50.619584 55.488272];
fi = [13.395837 12.373351 19.381854 26.244260 28.787526];

global xi yi zi;

xi = r.*cos(deg2rad(theta)).*cos(deg2rad(fi));
yi = r.*cos(deg2rad(theta)).*sin(deg2rad(fi));
zi = r.*sin(deg2rad(theta));

ti = [6.681975504357193e-02 6.684736756249902e-02 6.676997507374220e-02 6.675806218267925e-02 6.687535792178521e-02];
% ti = [6.6819755e-02 6.6847367e-02 6.6769975e-02 6.6758062e-02 6.6875357e-02];
global vsti

vsti = vs .* ti;
% ansx = 3.7199e6;
% ansy = 1.4392e6;
% ansz = 4.9773e6;

x0 = [0 0 0];

OPTIONS = optimset('Algorithm', 'levenberg-marquardt', 'Jacobian', 'on');
% OPTIONS.TolFun = 1e-16;
% OPTIONS.TolX = 1e-16;
% OPTIONS.MaxIter = 400;

[P, resnorm, residual, exitflag, output, lambda, Jacobian] = lsqnonlin(@FUN, x0, [], [], OPTIONS);

x = P(1);
y = P(2);
z = P(3);
rn = sqrt(x^2+y^2+z^2);
fin = rad2deg(atan(y/x));
thetan = rad2deg(asin(z/rn));

function [k, J] = FUN(x)
    global xi yi zi vsti;

    k = (x(1) - xi).^2 + (x(2) - yi).^2 + (x(3) - zi).^2 - vsti.^2;
    if nargout > 1
        J = [2*x(1)-2.*xi; 2*x(2)-2.*yi; 2*x(3)-2.*zi].';
    end
   
end