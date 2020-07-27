% FUNCTION: INITIAL ESTIMATION OF STATE VARIABLES
function f = fcn1(d, al, be, del, ph, g, zt1, zt2)
eps = 1e-3;
a1 = exp(g);
a2 = exp(d);
a3 = (ph + a2) / (1 + a2);
a4 = exp(a3 * g) / (al * be);
p = a4 - del / al;
a5 = a1 - del;
q = p - a5;
if (p < eps || q < eps)
    f = 1e12;
else
    lnp = log(p);
    lnq = log(q);
    lpq = lnp - lnq;
    a6 = log(1 - al);
    f1 = (zt1 - lpq)^2;
    f2 = (zt2 - d + a6)^2;
    f = f1 + f2;
end
