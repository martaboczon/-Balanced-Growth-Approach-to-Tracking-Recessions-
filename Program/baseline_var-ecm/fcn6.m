% FUNCTION: FORECASTING
function [f, crps_mat] =  fcn6(be, del, ph, T0, ny, npr, L, cp, s, zt, ...
    pi1, sigma1, pi2, sigma2)
ycn = cell(L, 1);
s = [s(T0 - 1 : T0, :); zeros(npr, size(s, 2))];
ztf = [zt(T0 - 1 : T0, :); zeros(npr, size(zt, 2))];
cpf = [cp(T0 - 1 : T0, :); zeros(npr, size(cp, 2))];
clf = zeros(npr, size(cp, 2)); 
zef = zeros(1, 7);
R1 = chol(sigma1, 'upper');
R2 = chol(sigma2, 'upper');
for l = 1 : L 
    for i = 1 : npr
        eps = randn(ny, 1);
        u = R1' * eps;
        j = i + 2;
        s(j,:) = pi1(1,:) + s(j-1,:) * pi1(2:4,:)...
            + s(j-2,:) * pi1(5:7,:) + u';
    end
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
    for i = 1 : npr
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
        eps = randn(ny, 1);
        u = R2' * eps;    
        ztf(j, 3 : 5) = zef * pi2 + u';
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
    end  
    ycn{l} = clf;
end
ycn = cell2mat(ycn);
ycn = reshape(ycn, npr, L, ny);
crps_mat = ycn;
av = mean(ycn, 2);
p25 = prctile(ycn, 2.5, 2);
p75 = prctile(ycn, 97.5, 2);
f = [av, p25, p75];
