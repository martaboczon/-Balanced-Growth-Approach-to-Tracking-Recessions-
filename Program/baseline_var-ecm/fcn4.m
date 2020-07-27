% FUNCTION: ERROR CORRECTION MECHANISM
function [X, Y] = fcn4(s, zt, T0, be, del, ph)
g = s(:, 1);
d = s(:, 2);
al = s(:, 3);
a1 = exp(g);
a2 = exp(d);
a3 = (ph + a2) ./ (1 + a2);
a4 = exp(a3 .* g) ./ (al * be);
p = a4 - del ./ al;
a5 = a1 - del;
q = p - a5;

time = 1 : T0;
t0 = time(5 : end);
t1 = time(4 : end - 1);
t2 = time(3 : end - 2);
Y = zeros(T0 - 4, 3);

% Bg's: k0
bg3k0 = al(t0)./(al(t0)-1).*log(p(t0))-al(t1)./(al(t1)-1)....
    .*log(p(t1));
bg4k0 = 1./(al(t0)-1).*log(p(t0))-1./(al(t1)-1).*log(p(t1))...
    +log(q(t0)./q(t1));
bg5k0 = d(t0)-d(t1)-log((1-al(t0))./(1-al(t1)))-log(p(t0)./p(t1))...
    +log(q(t0)./q(t1));

% Bg's: k1
bg3k1 = al(t1)./(al(t1)-1).*log(p(t1))-al(t2)./(al(t2)-1)...
    .*log(p(t2));
bg4k1 = 1./(al(t1)-1).*log(p(t1))-1./(al(t2)-1).*log(p(t2))...
    +log(q(t1)./q(t2));
bg5k1 = d(t1)-d(t2)- log(1-al(t1)) + log(1-al(t2))...
    - log(p(t1)) + log(p(t2)) + log(q(t1)) - log(q(t2));

% X Block:
X1 = ones(T0 - 4, 1);
X2 = zt(t1, 1) - log(p(t1) ./ q(t1));
X3 = zt(t1, 3) - g(t1) - bg3k1 ;
X4 = zt(t1, 4) - g(t1) - bg4k1;
X5 = zt(t1, 5) - bg5k1; 
X6 = g(t0) - g(t1);
X7 = al(t0) - al(t1);
X = [X1, X2, X3, X4, X5, X6, X7];

% Y Block:
Y(:, 1) = zt(t0, 3) - g(t0) - bg3k0 ;
Y(:, 2) = zt(t0, 4) - g(t0) - bg4k0;
Y(:, 3) = zt(t0, 5) - bg5k0;
