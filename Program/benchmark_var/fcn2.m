%Tracking
function clf =  fcn2(T0, npr, cp, zt, pi)
ztf = [zt(T0 - 2 : T0, :); zeros(npr, size(zt, 2))];
cpf = [cp(T0 - 2 : T0, :); zeros(npr, size(cp, 2))];
clf = zeros(1, size(cp, 2)); 
zef = zeros(1, 7);
        i = 1;
        j = i + 2;  
        zef(1) = 1;
        zef(2) = ztf(j - 1, 3) - ztf(j - 2, 3);
        zef(3) = ztf(j - 1, 4) - ztf(j - 2, 4);    
        zef(4) = ztf(j - 1, 5) - ztf(j - 2, 5);       
        zef(5) = ztf(j - 2, 3);     
        zef(6) = ztf(j - 2, 4);
        zef(7) = ztf(j - 2, 5);
        ztf(j, 3 : 5) = zef * pi;
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

