% Tracking
function clf =  fcn5(be, del, ph, T0, cp, s, zt, pi2)
s = s(T0 - 2 : T0, :);
ztf = zt(T0 - 2 : T0, :);
cpf = cp(T0 - 2 : T0, :);
clf = zeros(1, size(cp, 2)); 
zef = zeros(1, 7);
    g =  s(:, 1);
    d =  s(:, 2);
    al = s(:, 3);
    a1 = exp(g);
    a2 = exp(d);
    a3 = (ph + a2) ./ (1 + a2);
    a4 = exp(a3 .* g) ./ (al * be);
    p = a4 - del ./ al;
    a5 = a1 - del;
    q = p - a5;     
        i = 1;
        j = i + 2;     
        bgf3k0 = al(j)/(al(j)-1)*log(p(j))-al(j-1)/(al(j-1)-1)...
            *log(p(j-1));
        bgf4k0 = 1/(al(j)-1)*log(p(j))-1/(al(j-1)-1)*log(p(j-1))...
            +log(q(j)/q(j-1));
        bgf5k0 = d(j)-d(j-1)-log((1-al(j))/(1-al(j-1)))-log(p(j)/p(j-1))...
            +log(q(j)/q(j-1)); 
        bgf3k1 = al(j-1)/(al(j-1)-1)*log(p(j-1))-al(j-2)/(al(j-2)-1)...
            *log(p(j-2));
        bgf4k1 = 1/(al(j-1)-1)*log(p(j-1))-1/(al(j-2)-1)*log(p(j-2))...
            +log(q(j-1)/q(j-2));
        bgf5k1 = d(j-1)-d(j-2)-log((1-al(j-1))/(1-al(j-2)))...
            -log(p(j-1)/p(j-2))+log(q(j-1)/q(j-2));
        zef(1) = 1;
        zef(2) = ztf(j-1,1) - log(p(j-1) / q(j-1));
        zef(3) = ztf(j-1,3) - g(j-1) - bgf3k1;    
        zef(4) = ztf(j-1,4) - g(j-1) - bgf4k1;
        zef(5) = ztf(j-1,5) - bgf5k1;     
        zef(6) = s(j,1) - s(j-1,1);
        zef(7) = s(j,3) - s(j-1,3);
        ztf(j, 3 : 5) = zef * pi2;
        ztf(j, 3) = ztf(j, 3) + g(j) + bgf3k0;
        ztf(j, 4) = ztf(j, 4) + g(j) + bgf4k0;
        ztf(j, 5) = ztf(j, 5) + bgf5k0;
        ztf(j, 1) = ztf(j - 1, 1) + ztf(j, 3) - ztf(j, 4);
        ztf(j, 2) = ztf(j - 1, 2) + ztf(j, 3) - ztf(j, 4) + ztf(j, 5);
        cpf(j, :) = cpf(j - 1, :) + ztf(j, 3 : 5);
        clf(i, :) = exp(cpf(j, :));
        clf(i, 3) = 1 / (1 + clf(i, 3));
        clf(i, 1) = clf(i, 1) * clf(i, 3);
        clf(i, 2) = clf(i, 2) * clf(i, 3); 
