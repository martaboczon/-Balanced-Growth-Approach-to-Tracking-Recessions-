% Forecasting
function [f, crps_mat] =  fcn3(T0, ny, npr, L, cp, zt, ...
    pi, sigma)
ycn = cell(L, 1);
ztf = [zt(T0 - 1 : T0, :); zeros(npr, size(zt, 2))];
cpf = [cp(T0 - 1 : T0, :); zeros(npr, size(cp, 2))];
clf = zeros(npr, size(cp, 2)); 
zef = zeros(1, 7);
R2 = chol(sigma, 'upper');
for l = 1 : L 
    for i = 1 : npr
        j = i + 2;  
        zef(1) = 1;
        zef(2) = ztf(j - 1, 3) - ztf(j - 2, 3);
        zef(3) = ztf(j - 1, 4) - ztf(j - 2, 4);    
        zef(4) = ztf(j - 1, 5) - ztf(j - 2, 5);       
        zef(5) = ztf(j - 2, 3);     
        zef(6) = ztf(j - 2, 4);
        zef(7) = ztf(j - 2, 5);
        eps = randn(ny, 1);
        u = R2' * eps;    
        ztf(j, 3 : 5) = zef * pi + u';
        ztf(j, 3) = ztf(j, 3) + ztf(j - 1, 3);
        ztf(j, 4) = ztf(j, 4) + ztf(j - 1, 4);
        ztf(j, 5) = ztf(j, 5) + ztf(j - 1, 5);
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
