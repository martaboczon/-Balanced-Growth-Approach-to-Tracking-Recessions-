% FUNCTION: INITIAL ESTIMATION OF STATE VARIABLES
function f = fcn2(d1, al1, g1, st, be, del, ph, zt)
eps = 1e-3;
g = [st(1), g1];
d = [st(2), d1];
al = [st(3), al1];
a1 = exp(g);
a2 = exp(d);
a3 = (ph + a2) ./ (1 + a2);
a4 = exp(a3 .* g) ./ (al * be);
p = a4 - del ./ al;
a5 = a1 - del;
q = p - a5;
if (p(2) < eps || q(2) < eps)
    f = 1e12;
else
    lnp = log(p);
    lnq = log(q);
    a6 = log(1 - al);
    a7 = al ./ (al - 1); 
    a8 = 1 ./ (al - 1); 
    dlp = lnp(2) - lnp(1);
    dlq = lnq(2) - lnq(1);
    ztm(1) = lnp(2) - lnq(2);
    ztm(2) = d(2) - a6(2);
    ztm(3) = g(2) + a7(2) * lnp(2) - a7(1) * lnp(1);
    ztm(4) = g(2) + a8(2) * lnp(2) - a8(1) * lnp(1) + dlq;
    ztm(5) = (d(2) - d(1)) - (a6(2) - a6(1)) - dlp + dlq;
    f = sum((zt - ztm).^2);
end
