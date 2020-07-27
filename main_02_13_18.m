clc
clear
close all
rng(123)
%% Import data
datasource = 'fort331';
version = 4;
filename = sprintf('OutputF%i.xlsx', version);
dir = sprintf('C:\\Users\\marta\\Desktop\\RA\\SPR_2018\\EstimationF%i',...
    version);
cd(dir)
y = xlsread(datasource,'fort','B1:B1000');          %y
c = xlsread(datasource, 'fort', 'C1:C1000');        %c
n = xlsread(datasource, 'fort', 'F1:F1000');        %n
aFRED = xlsread(datasource, 'fort', 'G1:G1000');    %alpha
stopt = struct2cell(load('stopt.mat'));
stopt = cell2mat(stopt);
myrecession = [...
    datenum('Q4-1948', 'QQ-YY'), datenum('Q4-1949', 'QQ-YY');...
    datenum('Q2-1953', 'QQ-YY'), datenum('Q2-1954', 'QQ-YY');...
    datenum('Q3-1957', 'QQ-YY'), datenum('Q2-1958', 'QQ-YY');...
    datenum('Q2-1960', 'QQ-YY'), datenum('Q1-1961', 'QQ-YY');...
    datenum('Q4-1969', 'QQ-YY'), datenum('Q4-1970', 'QQ-YY');...
    datenum('Q4-1973', 'QQ-YY'), datenum('Q1-1975', 'QQ-YY');...
    datenum('Q1-1980', 'QQ-YY'), datenum('Q3-1980', 'QQ-YY');...
    datenum('Q3-1981', 'QQ-YY'), datenum('Q4-1982', 'QQ-YY');...
    datenum('Q3-1990', 'QQ-YY'), datenum('Q1-1991', 'QQ-YY');...
    datenum('Q1-2001', 'QQ-YY'), datenum('Q4-2001', 'QQ-YY');...
    datenum('Q4-2007', 'QQ-YY'), datenum('Q2-2009', 'QQ-YY')];
%% Initial parameters
warning('off','MATLAB:xlswrite:AddSheet')
be0 = 0.97;
del0 = 0.98;
ph0 = 1.3;
ts = 1.96;  % t-statistic threshold
npr = 4;    % forecasting horizon
L = 1000;   % number of generated forecasts
t1 = 230;
t2 = 270;
t = t1 : t2;
nt = size(t, 2);
out = cell(nt, 12);
fy = zeros(nt, 6, npr);
fc = zeros(nt, 6, npr);
fn = zeros(nt, 6, npr);
ft = zeros(nt, 3);
crps_mat = cell(nt, 1);
h = cell(nt, 1);
pValue = cell(nt, 1);
stat = cell(nt, 1);
cValue = cell(nt, 1);
mles = cell(nt, 1);
johansen = zeros(nt, 2);
%% Figure setting
startDate = datenum('Q1-1948', 'QQ-YYYY');
endDate = datenum('Q2-2019', 'QQ-YYYY');
xax = linspace(startDate, endDate, size(y, 1))';

startDate1 = datenum('Q2-1948', 'QQ-YYYY');
endDate1 = datenum('Q2-2019', 'QQ-YYYY');
xax1 = linspace(startDate1, endDate1, size(y, 1) - 1)';

sD0 = datenum('Q2-2005', 'QQ-YYYY');
eD0 = datenum('Q2-2015', 'QQ-YYYY');
fxax0 = linspace(sD0, eD0, nt)';

sD1 = datenum('Q3-1985', 'QQ-YYYY');
eD1 = datenum('Q3-2019', 'QQ-YYYY');
fxax1 = linspace(sD1, eD1, nt)';

sD2 = datenum('Q4-1985', 'QQ-YYYY');
eD2 = datenum('Q4-2019', 'QQ-YYYY');
fxax2 = linspace(sD2, eD2, nt)';

sD3 = datenum('Q1-1986', 'QQ-YYYY');
eD3 = datenum('Q1-2020', 'QQ-YYYY');
fxax3 = linspace(sD3, eD3, nt)';

sD4 = datenum('Q2-1986', 'QQ-YYYY');
eD4 = datenum('Q2-2020', 'QQ-YYYY');
fxax4 = linspace(sD4, eD4, nt)';
fxax = [fxax1, fxax2, fxax3, fxax4];
%% Recursive estimation
for i = 1 : nt
    pcav = pca([c(2 : t(i)) ./ c(1 : t(i) - 1),...
        y(2 : t(i)) ./ y(1 : t(i) - 1)]);
    c1 = pcav(1, 1);
    c2 = pcav(2, 1);
    au = c1 + c2;
    c1 = c1 / au;
    c2 = c2 / au;
    cl = [y(1 : t(i)) ./ n(1 : t(i)),...
          c(1 : t(i)) ./ n(1 : t(i)),...
         (1 - n(1 : t(i))) ./ n(1 : t(i))];
    cp = log(cl);
    zt1 = cp(:, 1) - cp(:, 2);
    zt2 = zt1 + cp(:, 3);
    zt3 = cp(2 : end, 1) - cp(1 : end - 1, 1);
    zt4 = cp(2 : end, 2) - cp(1 : end - 1, 2);
    zt5 = cp(2 : end, 3) - cp(1 : end - 1, 3);
    zt = [zt1(2 : end), zt2(2 : end), zt3, zt4, zt5];
    b0 = [zt3, zt4];
    st = [c1 * b0(:, 1) + c2 * b0(:, 2),...
        repmat(zeros(size(zt, 1), 1), 1, 2)];
    cp = cp(2 : end, :);
    T0 = size(st, 1);
    % Johansen cointegration test
%     [h{i},pValue{i},stat{i},cValue{i},mles{i}] =...
%         jcitest([y(1 : t(i)), c(1 : t(i))], 'model', 'H1*');
%     johansen(i, :) = table2array(h{i});
    % First period d and alpha
%     d0 = 1.40; 
%     al0 = 0.24;
%     x0 = [d0, al0];
%     OptSet = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,...
%         'Display', 'iter', 'TolFun', 1e-12,'TolX', 1e-12);
%     th = fminsearch(@(x) fcn1(x(1), x(2),...
%         be0, del0, ph0, st(1, 1), zt1(2), zt2(2)), x0, OptSet);
%     st(1, 2 : 3) = th;
%     % Subsequent periods d and alpha
%     d0 = 1.40; 
%     al0 = 0.24;
%     x0 = [d0, al0];
%     lb = [0.5, 0.2];
%     ub = [2, 0.4];
%     OptSet = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,...
%         'Display', 'iter', 'TolFun', 1e-12,'TolX', 1e-12);
%     for j = 2 : T0
%         th = fminsearch(@(x) fcn2(x(1), x(2),...
%             st(j, 1), st(j - 1, :), be0, del0, ph0, zt(j, :)), x0, OptSet);    
%         st(j, 2 : 3) = th; 
%     end  
    st = stopt(1 : T0, :);  
%     % Normal
%     if T0 == 241
%         st(241, 3) = st(241, 3) + 0.001;
%     end
%     if T0 == 242
%         st(242, 3) = st(242, 3) + 0.002;
%     end
%     if T0 == 243
%         st(243, 3) = st(243, 3) + 0.002;
%     end
%     if T0 == 244
%         st(244, 3) = st(244, 3) + 0.004;
%     end
%     if T0 == 245
%         st(245, 3) = st(245, 3) + 0.006;
%     end
%     if T0 == 246
%         st(246, 3) = st(246, 3) + 0.006;
%     end
%     if T0 == 247
%         st(247, 3) = st(247, 3) + 0.006;
%     end
%     if T0 == 248
%         st(248, 3) = st(248, 3) + 0.003; 
%     end
%     if T0 == 249
%         st(249, 3) = st(249, 3) + 0.003; 
%     end
%     % Still okay
%     if T0 == 245
%         st(245, 3) = st(245, 3) + 0.006; 
%     end
%     if T0 == 246
%         st(246, 3) = st(246, 3) + 0.006; 
%     end
%     if T0 == 247
%         st(247, 3) = st(247, 3) + 0.006;
%     end
%     if T0 == 248
%         st(248, 3) = st(248, 3) + 0.003; 
%     end
%     if T0 == 249
%         st(249, 3) = st(249, 3) + 0.003; 
%     end    
    %% Too late
    if T0 == 246
        st(246, 3) = st(246, 3) + 0.006; 
    end
    if T0 == 247
        st(247, 3) = st(247, 3) + 0.006;
    end
    if T0 == 248
        st(248, 3) = st(248, 3) + 0.003; 
    end
    if T0 == 249
        st(249, 3) = st(249, 3) + 0.003; 
    end        
    % Unrestricted SURE for (g, d, a) on constant c
    UNCON_sigma = cov(st(3 : end, :), 1);
    UNCON_pi = mean(st(3 : end, :))';
    UNCON_std = UNCON_pi ./ (sqrt(diag(UNCON_sigma)) / sqrt(T0 - 2));
    out{i, 1} = [UNCON_pi, UNCON_std];
    out{i, 2} = UNCON_sigma;
    out{i, 3} = [];   
    % Unrestricted SURE for (g, d, a) on (c, g-1, d-1, a-1, g-2, d-2, a-2)
    Y = st(3 : end, :);
    X = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)]; 
    ny = size(Y, 2);
    nx = size(X, 2);
    T = size(X, 1);    
    [USURE_pi, USURE_sigma, USURE_std] = fcn3(X, Y, nx, ny);
    out{i, 4} = [USURE_pi, USURE_std];
    out{i, 5} = USURE_sigma;
    out{i, 6} = [];       
    % Smoothed values
    X = blkdiag(X, X, X);
    nx = size(X, 2);
    s = X * reshape(USURE_pi, nx, 1);
    s = [zeros(2, 3); reshape(s, T, ny)];
    out{i, 14} = diag(eye(3)) - diag(USURE_sigma) ./ diag(UNCON_sigma);
    % Restricted SURE for (g, d, a) on (c, g-1, d-1, a-1, g-2, d-2, a-2)
%     Y = [st(3 : end, 1); st(3 : end, 2); st(3 : end, 3)];
%     Xg = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     Xd = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     Xa = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     X = blkdiag(Xg, Xd, Xa);
%     nx = size(Xg, 2) + size(Xd, 2) + size(Xa, 2);
%     T = size(Xg, 1);
%     [RSURE_pi, RSURE_sigma, RSURE_std, ~] = fcn4(X, Y, T, nx, ny, ts);
%     % Restricted SURE for (g, d, a) on (Dg-1, Dd-1, Da-1, g-2, d-2, a-2)
%     Y = [st(3 : end, 1); diff(st(2 : end, 2)); diff(st(2 : end, 3))];
%     Xg = [ones(T0 - 2, 1), diff(st(1 : end - 1, :)), st(1 : end - 2, :)];
%     Xd = [ones(T0 - 2, 1), diff(st(1 : end - 1, :)), st(1 : end - 2, :)];
%     Xa = [ones(T0 - 2, 1), diff(st(1 : end - 1, :)), st(1 : end - 2, :)];
%     X = blkdiag(Xg, Xd, Xa);
%     nx = size(X, 2);
%     T = size(Xg, 1);
%     [SURE1_pi, SURE1_sigma, SURE1_std, SURE1_xsx] = ...
%         fcn4(X, Y, T, nx, ny, ts);
%     [SURE2_pi, SURE2_std] = fcn5(SURE1_pi, SURE1_xsx, nx, ny);
%     Xg = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     Xd = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     Xa = [ones(T0 - 2, 1), st(2 : end - 1, :), st(1 : end - 2, :)];
%     out{i, 7} = [SURE1_pi, SURE1_std];
%     out{i, 8} = [SURE2_pi, SURE2_std];   
%     out{i, 9} = SURE1_sigma;
%     out{i, 10} = [];
    % Unconditional covariance matrix of dy = [zt3, zt4, zt5]
    dy_sigma = cov(zt(4 : end, 3 : 5), 1);
    dy_mean = mean(zt(4 : end, 3 : 5));
    % Unrestircted Error Correction Mechanism
    [X, Y] = fcn6(st, zt, T0, be0, del0, ph0);  
    ny = size(Y, 2);
    nx = size(X, 2);
    [ECM_pi, ECM_sigma, ECM_std] = fcn3(X, Y, nx, ny);
    out{i, 11} = [ECM_pi, ECM_std];
    out{i, 12} = ECM_sigma;      
    % Forecast
%     [f, crps_mat{i}] =  fcn7(be0, del0, ph0, T0, ny, npr, L, cp, st, zt,...
%         USURE_pi, USURE_sigma, ECM_pi, ECM_sigma);
    ft(i, :) = fcn11(be0, del0, ph0, T0, ny, npr, L, cp, st, zt,...
        USURE_pi, USURE_sigma, ECM_pi, ECM_sigma);
%     for j = 1 : npr     
%         fy(i, 1 : 5, j) = [i, T0 + j, f(j, :, 1)]; 
%         fc(i, 1 : 5, j) = [i, T0 + j, f(j, :, 2)]; 
%         fn(i, 1 : 5, j) = [i, T0 + j, f(j, :, 3)];              
%     end
end
dt = [y(t1 : t2), c(t1 : t2), n(t1 : t2)];
% for j = 1 : npr
%     fy(:,:,j) = [fy(:,1:2,j), [y(t1+j:end); zeros(j,1)], fy(:,3:5,j)]; 
%     fc(:,:,j) = [fc(:,1:2,j), [c(t1+j:end); zeros(j,1)], fc(:,3:5,j)]; 
%     fn(:,:,j) = [fn(:,1:2,j), [n(t1+j:end); zeros(j,1)], fn(:,3:5,j)]; 
% end
ft(ft==0)=NaN;
save('ft_toolate.mat', 'ft')
%save('ft.mat', 'ft')
% fy(fy==0)=NaN;
% fc(fc==0)=NaN;
% fn(fn==0)=NaN;
% %% Forecast accuracy
% load('C:\Users\marta\Desktop\RA\SPR_2018\EstimationF3\var.mat')
% cd 'C:\Users\marta\Desktop\RA\SPR_2018\EstimationF1'
% nber(:, :, 1) = xlsread('nberdates','nber1','E2:AZ1000');   
% nber(:, :, 2) = xlsread('nberdates','nber2','E2:AZ1000');  
% nber(:, :, 3) = xlsread('nberdates','nber3','E2:AZ1000');  
% nber(:, :, 4) = xlsread('nberdates','nber4','E2:AZ1000');
% nber(:, :, 5) = xlsread('nberdates','nber0','E2:AZ1000');
% nber(nber==0)=NaN;
% for j = 1 : npr
%     xlswrite(filename, fy(:, :, j), ['fy' num2str(j)], 'A2')
%     xlswrite(filename, fc(:, :, j), ['fc' num2str(j)], 'A2')
%     xlswrite(filename, fn(:, :, j), ['fn' num2str(j)], 'A2')  
% end
% % Root mean square error (RMSE)
% fcn9(filename, npr, fy, fc, fn, fy_var, fc_var, fn_var, nber,...
%     dt, ft, ft_var)
% % Mean absolute error (MAE)
% fcn10(filename, npr, fy, fc, fn, fy_var, fc_var, fn_var, nber,...
%     dt, ft, ft_var)
% %Continous rank probability score (CRPS)
% [crps_ecm1, crps_ecm2, crps_ecm3, crps_var1, crps_var2, crps_var3] =... 
%     fcn12(filename, npr, nt, ny, L, fy, fc, fn, crps_mat, crps_var, nber);
% %% Output display
% % for i = 1 : T0 - t1 + 1 
% %     xlswrite(filename, [out{i, 1}, out{i, 2}], num2str(t1 + i - 1), 'A2')
% %     xlswrite(filename, [out{i, 4}], num2str(t1 + i - 1), 'A6')
% %     xlswrite(filename, [out{i, 5}], num2str(t1 + i - 1), 'A14')
% %     xlswrite(filename, [out{i, 7}], num2str(t1 + i - 1), 'A18')
% %     xlswrite(filename, [out{i, 8}], num2str(t1 + i - 1), 'A26')
% %     xlswrite(filename, [out{i, 9}], num2str(t1 + i - 1), 'A34')
% %     xlswrite(filename, [out{i, 11}], num2str(t1 + i - 1), 'A38')
% %     xlswrite(filename, [out{i, 12}], num2str(t1 + i - 1), 'A47')
% % end
% %% Balanced growth ratios
% fig = figure;
% left_color = [0 0 0];
% right_color = [0 0 1];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% yyaxis left
%     plot(xax1, zt(:, 1), 'k:', 'LineWidth', 1);   
%     xlabel('Year')
%     set(gca,'yLim',[0.1, 0.4],'YTick', 0.1 : 0.1 : 0.4)      
% yyaxis right
%     plot(xax1, zt(:, 2), 'b', 'LineWidth', 1);   
%     set(gca,'yLim',[0.6, 0.9],'YTick', 0.6 : 0.1 : 0.9)  
%     xlabel('Year')
%     datetick('x','YYYY');
%     set(gca, 'FontSize', 12)
%     set(gca, 'FontName', 'Times New Roman')
%     hBand = recessionplot('recessions', myrecession);   
%     set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
%     leg = legend({'$$\log (y_t / c_t)$$',...
%         '$$\log [y_t/c_t \times (1-n_t) /n_t] $$'}, 'interpreter', 'latex');
%     set(leg, 'Location', 'northwest')
%     legend boxoff
% saveas(gcf, 'series_BG2.png')
% %% Laws of motion
% figure
%     subplot(3, 1, 1)
%         plot(xax1, zt(:, 3), 'b', 'LineWidth', 1);  
%         hline(0, 'k:')        
%         xlabel('Year')
%         set(gca,'yLim',[-0.03, 0.05],'YTick', -0.03 : 0.02 : 0.05)      
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12)
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot;   
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)       
%         xlabel('Year')
%         title('$$\Delta \log (y_t / n_t)$$', 'interpreter', 'latex')
%     subplot(3, 1, 2)
%         plot(xax1, zt(:, 4), 'b', 'LineWidth', 1); 
%         hline(0, 'k:')        
%         xlabel('Year')
%         set(gca,'yLim',[-0.03, 0.05],'YTick', -0.03 : 0.02 : 0.05)      
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12)
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot;   
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)        
%         xlabel('Year')
%         title('$$\Delta \log (c_t / n_t)$$', 'interpreter', 'latex')        
%     subplot(3, 1, 3)
%         plot(xax1, zt(:, 5), 'b', 'LineWidth', 1);
%         hline(0, 'k:')
%         xlabel('Year')
%         set(gca,'yLim', [-0.05, 0.07],'YTick', -0.05 : 0.03 : 0.07)      
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12)
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot;   
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)       
%         xlabel('Year')
%         title('$$\Delta \log (1-n_t / n_t)$$', 'interpreter', 'latex')        
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperType', 'A4');
%     size = get(gcf, 'PaperSize');
%     width = size(1);
%     height = size(2);
%     set(gcf, 'PaperPosition', [0, 0, width, height])
%     print(gcf,'laws_of_motion','-dpdf')          
% %% y, c, and n series
% figure
% 	subplot(3, 1, 1)
%         plot(xax, y / 10000, 'b', 'LineWidth', 1.5);
%         axis([-Inf Inf 0.5 3.5])
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12) 
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot; 
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
%         set(gca,'yLim',[0.5, 3.5],'YTick', 0.5 : 1 : 3.5)     
%         xlabel('Year')
%         ylabel('10,000 of chained 2012 dollars')  
%         title('\it Real output per capita', 'FontWeight','Normal')
%     subplot(3, 1, 2)	
%         plot(xax, c / 10000, 'b', 'LineWidth', 1.5);
%         axis([-Inf Inf 0.5 3.5])    
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12) 
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot;
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
%         set(gca,'yLim',[0.5, 3.5],'YTick', 0.5 : 1 : 3.5)      
%         xlabel('Year')
%         ylabel('10,000 of chained 2012 dollars')
%         title('\it Real consmuption per capita', 'FontWeight','Normal')     
%     subplot(3, 1, 3)
%         plot(xax, n * 100, 'b', 'LineWidth', 1.5);
%         axis([-Inf Inf 32 41])
%         datetick('x','YYYY');
%         set(gca, 'FontSize', 12) 
%         set(gca, 'FontName', 'Times New Roman')
%         hBand = recessionplot;
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
%         xlabel('Year') 
%         ylabel('Percentage')
%         title('\it Fraction of time spent working', 'FontWeight','Normal')        
%         set(gca,'yLim',[32, 41],'YTick', 32 : 3 : 41)  
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperType', 'A4');
%     size = get(gcf, 'PaperSize');
%     width = size(1);
%     height = size(2);
%     set(gcf, 'PaperPosition', [0, 0, width, height])
%     print(gcf,'ycn','-dpdf')        
% %% State variable series
% figure
%     plot(xax1, st(:, 1), 'k:', 'LineWidth', 0.5);
%     hold on
%     plot(xax1(3 : end), s(3 : end, 1), 'b', 'LineWidth', 1.5);
%     hline(0, 'k:')
%     axis([-Inf Inf -0.035 0.045])
%     ax = gca;
%     ax.XTick = xax;
%     datetick('x','YYYY');
%     set(gca, 'FontName', 'Times New Roman')
%     hBand = recessionplot;
%     set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
%     set(gca, 'FontSize', 14) 
%     xlabel('Year')
%     set(gca, 'yLim', [-0.02, 0.03],'YTick', -0.02 : 0.01 : 0.03)       
%     leg = legend({'Initial values', 'Fitted values'},...
%         'interpreter', 'latex');
%     set(leg, 'Location', 'southwest', 'FontSize', 12)
%     legend boxoff
%     saveas(gcf, 'series_G.png')     
% figure
%     plot(xax1, 1 ./ (exp(st(:, 2)) + 1), 'k:', 'LineWidth', 0.5);
%     hold on
%     plot(xax1(3 : end), 1 ./ (exp(s(3 : end, 2)) + 1) ,...
%         'b', 'LineWidth', 1.5); 
%     axis([-Inf Inf 0.37 0.47])    
%     ax = gca;
%     ax.XTick = xax;
%     datetick('x','YYYY');
%     set(gca, 'FontName', 'Times New Roman')
%     hBand = recessionplot('recessions', myrecession);  
%     set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)        
%     xlabel('Year')     
%     set(gca, 'FontSize', 14) 
%     xlabel('Year')
%     set(gca, 'yLim', [0.37, 0.47],'YTick', 0.37 : 0.02 : 0.47)  
%     leg = legend({'Initial values', 'Fitted values'},...
%         'interpreter', 'latex');
%     set(leg, 'Location', 'southwest', 'FontSize', 12)
%     legend boxoff    
%     saveas(gcf, 'series_D.png')      
% figure
%     plot(xax1, st(:, 3), 'k:', 'LineWidth', 0.5);
%     hold on
%     plot(xax1(3 : end), s(3 : end, 3), 'b', 'LineWidth', 1.5);  
%     hold on
%     plot(xax, aFRED, 'r-', 'LineWidth', 1);
%     axis([-Inf Inf 0.26 0.46])    
%     ax = gca;
%     ax.XTick = xax;
%     datetick('x','YYYY');
%     set(gca, 'yLim', [0.26, 0.46],'YTick', 0.26 : 0.04 : 0.46)      
%     set(gca, 'FontSize', 14) 
%     set(gca, 'FontName', 'Times New Roman')
%     hBand = recessionplot('recessions', myrecession); 
%     set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
%     xlabel('Year')      
%     leg = legend({'Initial values', 'Fitted values', 'Data'},...
%         'interpreter', 'latex');
%     set(leg, 'Location', 'northwest', 'FontSize', 12)
%     legend boxoff    
%     saveas(gcf, 'series_A.png')   
% %% VAR parameter stability
% b = cell2mat(out(:, 4));
% ax1 = 1.96;
% figure
% k = 0;  
% for ii = 1 : 7
%     for j = 1 : 3   
%         k = k + 1;
%         subplot(7, 3, k)
%         bj = b(ii : 7 : end, j);
%         tj = b(ii : 7 : end, j + 3);
%         sj = bj ./ tj;
%         x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
%         inBetween = [(bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
%         fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
%         hold on
%         plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)
%         hold on 
%         plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
%         hold on 
%         plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5) 
%         hline(0, 'k:')
%         ax = gca;
%         ax.XTick = xax;
%         if ii == 1
%             set(gca,'yLim',[-0.3, 0.3],'YTick', -0.3 : 0.3 : 0.3)
%         elseif ii == 2
%             set(gca,'yLim',[0, 10],'YTick', 0 : 5 : 10)          
%         elseif ii == 3
%             set(gca,'yLim',[-0.5, 1.5],'YTick', -0.5 : 1 : 1.5)    
%         elseif ii == 4
%             set(gca,'yLim',[-4, 0],'YTick', -4 : 2 : 0)     
%         elseif ii == 5
%             set(gca,'yLim',[-10, 0],'YTick', -10 : 5 : 0)    
%         elseif ii == 6
%             set(gca,'yLim',[-0.6, 0],'YTick', -0.6 : 0.3 : 0)  
%         elseif ii == 7
%             set(gca,'yLim',[0, 5],'YTick', 0 : 2.5 : 5)                
%         end            
%         datetick('x','YYYY');  
%         ax.XTickLabel = {'1980', '', '2000', '', '2020'};         
%         set(gca, 'FontSize', 10)          
%         hBand = recessionplot('recessions', myrecession); 
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)
%         set(gca, 'FontName', 'Times New Roman')
%     end       
% end    
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperType', 'A4');
% size = get(gcf, 'PaperSize');
% width = size(1);
% height = size(2);
% set(gcf, 'PaperPosition', [0, 0, width, height])
% print(gcf,'var','-dpdf')  
% %% ECM parameter stability
% b = cell2mat(out(:, 11));
% ax1 = 1.96;
% figure
% k = 0;   
% for ii = 1 : 7
%     for j = 1 : 3  
%         k = k + 1;
%         subplot(7, 3, k)
%         bj = b(ii : 7 : end, j);
%         tj = b(ii : 7 : end, j + 3);
%         sj = bj ./ tj;
%         x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
%         inBetween = [(bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
%         fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
%         hold on
%         plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)
%         hold on 
%         plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
%         hold on 
%         plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5)    
%         hline(0, 'k:')
%         ax = gca;
%         ax.XTick = xax;
%         if ii == 1
%             set(gca,'yLim',[-0.003, 0.003],'YTick', -0.003 : 0.003 : 0.003)
%         elseif ii == 2
%             set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.15 : 0.15)          
%         elseif ii == 3
%             set(gca,'yLim',[-1, 2],'YTick', -1 : 1.5 : 2)    
%         elseif ii == 4
%             set(gca,'yLim',[-1, 1.5],'YTick', -1 : 1.25 : 1.5)     
%         elseif ii == 5
%             set(gca,'yLim',[-0.2, 0.8],'YTick', -0.2 : 0.5 : 0.8)    
%         elseif ii == 6
%             set(gca,'yLim',[-1, 6],'YTick', -1 : 3.5 : 6)  
%         elseif ii == 7
%             set(gca,'yLim',[-2, 2],'YTick', -2 : 2 : 2)                
%         end            
%         datetick('x','YYYY');  
%         ax.XTickLabel = {'1980', '', '2000', '', '2020'};         
%         set(gca, 'FontSize', 10)          
%         hBand = recessionplot('recessions', myrecession);
%         set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)
%         set(gca, 'FontName', 'Times New Roman')
%     end      
% end
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperType', 'A4');
% size = get(gcf, 'PaperSize');
% width = size(1);
% height = size(2);
% set(gcf, 'PaperPosition', [0, 0, width, height])
% print(gcf,'ecm','-dpdf') 
% %% ECM components
% b = cell2mat(out(:, 11));
% fig = figure;
% for j = 1 : 3     
%     subplot(3, 1, j)
%         bj = b(2 : 7 : end, j);
%         tj = b(2 : 7 : end, j + 3);
%         sj = bj ./ tj;
%         x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
%         inBetween = [ (bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
%         fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
%         hold on
%         plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)  
%         hold on
%         plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
%         hold on 
%         plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5)            
%         hline(0, 'k:')
%         ax = gca;
%         ax.XTick = xax;
%         yt = yticks;      
%         datetick('x','YYYY');  
%         set(gca, 'FontSize', 14)
%         if j == 1
%             set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)  
%         elseif j == 2
%             set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)   
%         elseif j == 3
%             set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)              
%         end
%         set(gca, 'FontName', 'Times New Roman')        
%         hBand = recessionplot('recessions', myrecession); 
%         set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)  
%         xlabel('Year')       
% end  
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperType', 'A4');
% size = get(gcf, 'PaperSize');
% width = size(1);
% height = size(2);
% set(gcf, 'PaperPosition', [0, 0, width, height])
% print(gcf,'ecmcom','-dpdf')
%% Figures 1-4
ft = struct2cell(load('ft.mat'));
ft = cell2mat(ft);
ft_normal = struct2cell(load('ft_normal.mat'));
ft_normal = cell2mat(ft_normal);

ft_stillokay = struct2cell(load('ft_stillokay.mat'));
ft_stillokay = cell2mat(ft_stillokay);

ft_toolate = struct2cell(load('ft_toolate.mat'));
ft_toolate = cell2mat(ft_toolate);
figure
j = 0;
for k = 1 : 3
        if k == 1
            ftshock = ft_normal;
        elseif k == 2
            ftshock = ft_stillokay;
        elseif k == 3
            ftshock = ft_toolate;
        end
    for i = 1 : 2
        j = j + 1;
        if i == 1          
            z = fy / 10000;
    %         zvar = fy_var / 10000;
            v = 'Y';
            ftt = ft(:, 1) / 10000;
            fttshock = ftshock(:, 1) / 10000;
            dtt = dt(:, 1) / 10000;
    %         ztt = ft_var(:, 1) / 10000;        
        elseif i == 2
            z = fc / 10000;
    %         zvar = fc_var / 10000;        
            v = 'C';  
            ftt = ft(:, 2) / 10000;
            fttshock = ftshock(:, 2) / 10000;
            dtt = dt(:, 2) / 10000;
    %         ztt = ft_var(:, 2) / 10000;        
        elseif i == 3
            z = fn * 100;
    %         zvar = fn_var * 100;        
            v = 'N';          
            ftt = ft(:, 3) * 100;
            fttshock = ftshock(:, 3) * 100;
            dtt = dt(:, 3) * 100;
    %         ztt = ft_var(:, 3) * 100;        
        end
    sD1 = datenum('01-01-1985');
    eD1 = datenum('10-01-2020');
    xaxl = linspace(sD1, eD1, 144)';
        subplot(3, 2, j)
        h1 = plot(fxax0, dtt);
        hold on
        h3 = plot(fxax0, ftt);  
        hold on
        h2 = plot(fxax0, fttshock);          
        hold on
        plot(xaxl, -100, 'w')
        set(h1, 'LineWidth', 2, 'Color', 'b')
        set(h3, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-.')
        set(h2, 'LineWidth', 1, 'Color', 'r')       
        if i == 1
            set(gca,'yLim',[2.2 2.8])  
            ylabel('10,000 of chained 2012 dollars')
        elseif i == 2 
            set(gca,'yLim',[1.6 2.2])
            ylabel('10,000 of chained 2012 dollars')            
        elseif i == 3
            set(gca,'yLim',[32 38]) 
            ylabel('Percentage')            
        end    
        if i == 1 && k == 1
            title('{\it Policy 1: Output (y)}', 'FontWeight','Normal')
        end
        if i == 2 && k == 1
            title('{\it Policy 1: Consumption (c)}', 'FontWeight','Normal')
        end
        if i == 1 && k == 2
            title('{\it Policy 2: Output (y)}', 'FontWeight','Normal')
        end
        if i == 2 && k == 2
            title('{\it Policy 2: Consumption (c)}', 'FontWeight','Normal')
        end  
        if i == 1 && k == 3
            title('{\it Policy 3: Output (y)}', 'FontWeight','Normal')
        end
        if i == 2 && k == 3
            title('{\it Policy 3: Consumption (c)}', 'FontWeight','Normal')
        end             
        ax = gca;
        ax.XTick = xax;
        ax.TitleFontSizeMultiplier = 1.5;        
        datetick('x','YYYY')
        xlim([datenum('04-01-2005') datenum('04-01-2015')])
        xticks([datenum('04-01-2005'), ...
            datenum('04-01-2007'),...
            datenum('04-01-2009'),...
            datenum('04-01-2011'),...
            datenum('04-01-2013'), ...
            datenum('04-01-2015')])
        xticklabels({'Q2-2005','Q2-2007','Q2-2009',...
            'Q2-2011','Q2-2013','Q2-2015'})
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot('recessions', myrecession);
        set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
        xlabel('Year')     
        leg = legend([h1, h2, h3],...
            {'Data', 'VAR-ECM fitted value: Policy experiment',...
            'VAR-ECM fitted value: Baseline model'});
        legend boxoff
        set(leg, 'Location', 'Northwest', 'FontSize', 7)
        set(gca, 'FontName', 'Times New Roman')  
        set(gca, 'FontSize', 9)        
     end      
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf,'poliexp','-dpdf')

% %% Figures recessions' windows
% close all
% load('C:\Users\marta\Desktop\RA\SPR_2018\EstimationF3\var.mat')
% cd 'C:\Users\marta\Desktop\RA\SPR_2018\EstimationF1'
% rec1(:, :, 1) = xlsread('recession1','nber1','E2:Z1000');   
% rec1(:, :, 2) = xlsread('recession1','nber2','E2:Z1000');  
% rec1(:, :, 3) = xlsread('recession1','nber3','E2:Z1000');  
% rec1(:, :, 4) = xlsread('recession1','nber4','E2:Z1000');
% rec1(:, :, 5) = xlsread('recession1','nber0','E2:Z1000');  
% 
% rec2(:, :, 1) = xlsread('recession2','nber1','E2:Z1000');   
% rec2(:, :, 2) = xlsread('recession2','nber2','E2:Z1000');  
% rec2(:, :, 3) = xlsread('recession2','nber3','E2:Z1000');  
% rec2(:, :, 4) = xlsread('recession2','nber4','E2:Z1000');
% rec2(:, :, 5) = xlsread('recession2','nber0','E2:Z1000');
% 
% rec3(:, :, 1) = xlsread('recession3','nber1','E2:Z1000');   
% rec3(:, :, 2) = xlsread('recession3','nber2','E2:Z1000');  
% rec3(:, :, 3) = xlsread('recession3','nber3','E2:Z1000');  
% rec3(:, :, 4) = xlsread('recession3','nber4','E2:Z1000');
% rec3(:, :, 5) = xlsread('recession3','nber0','E2:Z1000');
% 
% for i = 1 : 3
%     if i == 1
%         z = fy / 10000;
%         zvar = fy_var / 10000;
%         v = 'Y';
%         ftt = ft(:, 1) / 10000;
%         dtt = dt(:, 1) / 10000;
%         ztt = ft_var(:, 1) / 10000;
%     elseif i == 2
%         z = fc / 10000;
%         zvar = fc_var / 10000;        
%         v = 'C';        
%         ftt = ft(:, 2) / 10000;
%         dtt = dt(:, 2) / 10000;  
%         ztt = ft_var(:, 2) / 10000;        
%     else
%         z = fn * 100;
%         zvar = fn_var * 100;        
%         v = 'N';      
%         ftt = ft(:, 3) * 100;
%         dtt = dt(:, 3) * 100;  
%         ztt = ft_var(:, 3) * 100;        
%     end    
%     sD1 = datenum('01-01-1985');
%     eD1 = datenum('10-01-2020');
%     xaxl = linspace(sD1, eD1, 144)';
%     for r = 1 : 3
%         fig = figure('units','normalized','outerposition',[0 0 1 1]);       
%         if r == 1
%             rec = rec1;
%         elseif r == 2 
%             rec = rec2;
%         elseif r == 3
%             rec = rec3;
%         end
%         rec(rec==0)=NaN;
%         j = 0;
%         for s = 1 : 4
%             subplot(2, 2, s)
%             if s == 1
%                 h2 = plot(fxax0, ftt .* rec(:, 9, 5)); 
%                 hold on
%                 h3 = plot(fxax0, ztt .* rec(:, 9, 5)); 
%                 hold on
%                 h1 = plot(fxax0, dtt .* rec(:, 9, 5)  , '-');
%                 set(h1, 'LineWidth', 2.5, 'Color', 'b')
%                 set(h2, 'LineWidth', 1, 'Color', 'k')
%                 set(h3, 'LineWidth', 1, 'Color', [0.7 0.7 0.7])                
%             else      
%                 j = j + 1;
%                 h2 = plot(fxax(:, j), z(:, 4, j) .* rec(:, 9, j));  
%                 hold on
%                 h4 = plot(fxax(:, j), z(:, 5, j) .* rec(:, 9, j),...
%                     '--', 'Color', 'k', 'LineWidth', 0.5);  
%                 hold on
%                 plot(fxax(:, j), z(:, 6, j) .* rec(:, 9, j),...
%                     '--', 'Color', 'k', 'LineWidth', 0.5);  
%                 hold on
%                 h3 = plot(fxax(:, j), zvar(:, 4, j) .* rec(:, 9, j));  
%                 hold on
%                 h5 = plot(fxax(:, j), zvar(:, 5, j) .* rec(:, 9, j),...
%                     'Color', [0.5 0.5 0.5], 'LineStyle', '-.',...
%                     'LineWidth', 0.1);   
%                 hold on
%                 plot(fxax(:, j), zvar(:, 6, j) .* rec(:, 9, j),...
%                     'Color', [0.5 0.5 0.5], 'LineStyle', '-.',...
%                     'LineWidth', 0.1);    
%                 hold on
%                 h1 = plot(fxax(:, j), z(:, 3, j) .* rec(:, 9, j)  , '-');
%                 set(h1, 'LineWidth', 2.5, 'Color', 'b')
%                 set(h2, 'LineWidth', 1, 'Color', 'k')
%                 set(h3, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]) 
%             end
%             hold on
%             plot(xaxl, -100, 'w')
%             hold on
%             if i == 1 && r == 1
%                 set(gca,'yLim', [1.5, 1.8])  
%                 ylabel('10,000 of chained 2012 dollars')
%             elseif i == 1 && r == 2
%                 set(gca,'yLim', [1.9, 2.35])  
%                 ylabel('10,000 of chained 2012 dollars')
%             elseif i == 1 && r == 3
%                 set(gca,'yLim', [2.1, 2.6])  
%                 ylabel('10,000 of chained 2012 dollars')                  
%             elseif i == 2 && r == 1               
%                 set(gca,'yLim', [1.3, 1.5])
%                 ylabel('10,000 of chained 2012 dollars')      
%             elseif i == 2 && r == 2               
%                 set(gca,'yLim', [1.55, 1.85])
%                 ylabel('10,000 of chained 2012 dollars')
%             elseif i == 2 && r == 3               
%                 set(gca,'yLim', [1.8, 2.05])
%                 ylabel('10,000 of chained 2012 dollars')                   
%             elseif i == 3 && r == 1               
%                 set(gca,'yLim', [32, 38])
%                 ylabel('Percentage') 
%             elseif i == 3 && r == 2               
%                 set(gca,'yLim', [33, 39])
%                 ylabel('Percentage')            
%             elseif i == 3 && r == 3               
%                 set(gca,'yLim', [31, 37])
%                 ylabel('Percentage')                   
%             end                     
%             hBand = recessionplot('recessions', myrecession);
%             set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
%             ax = gca;
%             ax.XTick = xax;
%             ax.TitleFontSizeMultiplier = 1.5;        
%             datetick('x','QQ-YYYY')
%             if r == 1 
%                 xlim([datenum('Q1-1990', 'QQ-YY')...
%                       datenum('Q3-1992', 'QQ-YY')])
%                 xticks([datenum('Q1-1990', 'QQ-YY')...
%                         datenum('Q3-1990', 'QQ-YY')...
%                         datenum('Q1-1991', 'QQ-YY')...
%                         datenum('Q3-1991', 'QQ-YY')...
%                         datenum('Q1-1992', 'QQ-YY')...
%                         datenum('Q3-1992', 'QQ-YY')])
%                 xticklabels({'Q1-1990', 'Q3-1990', 'Q1-1991', 'Q3-1991',...
%                              'Q1-1992', 'Q3-1992'})                
%             elseif r == 2 
%                 xlim([datenum('Q3-2000', 'QQ-YY')...
%                       datenum('Q2-2003', 'QQ-YY')])
%                 xticks([datenum('Q3-2000', 'QQ-YY')...
%                       datenum('Q1-2001', 'QQ-YY')...
%                       datenum('Q3-2001', 'QQ-YY')...
%                       datenum('Q1-2002', 'QQ-YY')...
%                       datenum('Q3-2002', 'QQ-YY')...
%                       datenum('Q1-2003', 'QQ-YY')])
%                 xticklabels({'Q3-2000', 'Q1-2001', 'Q3-2001', 'Q1-2002',...
%                              'Q3-2002', 'Q1-2003'})                        
%             elseif r == 3 
%                 xlim([datenum('Q2-2007', 'QQ-YY')...
%                       datenum('Q4-2010', 'QQ-YY')])                
%                 xticks([datenum('Q2-2007', 'QQ-YY')...
%                       datenum('Q4-2007', 'QQ-YY')...
%                       datenum('Q2-2008', 'QQ-YY')...
%                       datenum('Q4-2008', 'QQ-YY')...
%                       datenum('Q2-2009', 'QQ-YY')...
%                       datenum('Q4-2009', 'QQ-YY')...
%                       datenum('Q2-2010', 'QQ-YY')...
%                       datenum('Q4-2010', 'QQ-YY')])
%                 xticklabels({'Q2-2007', 'Q4-2007', 'Q2-2008', 'Q4-2008',...
%                          'Q2-2009', 'Q4-2009', 'Q2-2010', 'Q4-2010'})                 
%             end
%             xlabel('Year')     
%             set(gca, 'FontName', 'Times New Roman')  
%             set(gca, 'FontSize', 10)   
%             if s == 1
%                 leg = legend([h1, h2, h3],...
%                      {'Data', 'VAR-ECM fitted value',...
%                      'Baseline VAR fitted value'});           
%                 title('\it 0-step ahead', 'FontWeight','Normal')                 
%             else
%                 leg = legend([h1, h2, h4, h3, h5],...
%                      {'Data', 'VAR-ECM point forecast',...
%                      '95% confidence interval for VAR-ECM forecast',...
%                      'Baseline VAR point forecast',...
%                      '95% confidence interval for baseline VAR forecast'});
%                 title(['\it' num2str(j) '-step ahead'],...
%                     'FontWeight', 'Normal')                 
%             end
%             legend boxoff
%             if i == 1 || i == 3
%                 set(leg, 'Location', 'Southwest', 'FontSize', 6)
%             else
%                 set(leg, 'Location', 'Northwest', 'FontSize', 6)
%             end
%             saveas(gcf, ['rec_' num2str(r) 'var_' num2str(i) '.png'])
%         end
%     end
% end
% 
% save('stopt.mat', 'st');
