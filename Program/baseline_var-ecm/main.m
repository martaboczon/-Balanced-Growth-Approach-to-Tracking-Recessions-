clc
clear
close all
rng(123)
%% Import data
datasource = [pwd '\input\data.xlsx'];
filename = [pwd '\output\output.xls'];
y = xlsread(datasource,'data','D2:D1000');          %y
c = xlsread(datasource, 'data', 'E2:E1000');        %c
n = xlsread(datasource, 'data', 'F2:F1000');        %n
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
t1 = 150;
t2 = size(y, 1);
t = t1 : t2;
nt = size(t, 2);
out = cell(nt, 12);
fy = zeros(nt, 6, npr);
fc = zeros(nt, 6, npr);
fn = zeros(nt, 6, npr);
ft = zeros(nt, 3);
crps_mat = cell(nt, 1);
%% Figure setting
startDate = datenum('Q1-1948', 'QQ-YYYY');
endDate = datenum('Q2-2019', 'QQ-YYYY');
xax = linspace(startDate, endDate, size(y, 1))';

startDate1 = datenum('Q2-1948', 'QQ-YYYY');
endDate1 = datenum('Q2-2019', 'QQ-YYYY');
xax1 = linspace(startDate1, endDate1, size(y, 1) - 1)';

sD0 = datenum('Q2-1985', 'QQ-YYYY');
eD0 = datenum('Q2-2019', 'QQ-YYYY');
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
    % First period d and alpha
    d0 = 1.40; 
    al0 = 0.24;
    x0 = [d0, al0];
    OptSet = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,...
        'Display', 'iter', 'TolFun', 1e-12,'TolX', 1e-12);
    th = fminsearch(@(x) fcn1(x(1), x(2),...
        be0, del0, ph0, st(1, 1), zt1(2), zt2(2)), x0, OptSet);
    st(1, 2 : 3) = th;
    % Subsequent periods d and alpha
    d0 = 1.40; 
    al0 = 0.24;
    x0 = [d0, al0];
    lb = [0.5, 0.2];
    ub = [2, 0.4];
    OptSet = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,...
        'Display', 'iter', 'TolFun', 1e-12,'TolX', 1e-12);
    for j = 2 : T0
        th = fminsearch(@(x) fcn2(x(1), x(2),...
            st(j, 1), st(j - 1, :), be0, del0, ph0, zt(j, :)), x0, OptSet);    
        st(j, 2 : 3) = th; 
    end  
    % VAR process
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
    % Error Correction Mechanism
    [X, Y] = fcn4(st, zt, T0, be0, del0, ph0);  
    ny = size(Y, 2);
    nx = size(X, 2);
    [ECM_pi, ECM_sigma, ECM_std] = fcn3(X, Y, nx, ny);
    out{i, 11} = [ECM_pi, ECM_std];
    out{i, 12} = ECM_sigma;    
    % Tracking
    ft(i, :) = fcn5(be0, del0, ph0, T0, cp, st, zt, ECM_pi);    
    % Forecasting
    [f, crps_mat{i}] =  fcn6(be0, del0, ph0, T0, ny, npr, L, cp, st, zt,...
        USURE_pi, USURE_sigma, ECM_pi, ECM_sigma);
    for j = 1 : npr     
        fy(i, 1 : 5, j) = [i, T0 + 1 + j, f(j, :, 1)]; 
        fc(i, 1 : 5, j) = [i, T0 + 1 + j, f(j, :, 2)]; 
        fn(i, 1 : 5, j) = [i, T0 + 1 + j, f(j, :, 3)];              
    end
end
dt = [y(t1 : end), c(t1 : end), n(t1 : end)];
for j = 1 : npr
    fy(:,:,j) = [fy(:,1:2,j), [y(t1+j:end); zeros(j,1)], fy(:,3:5,j)]; 
    fc(:,:,j) = [fc(:,1:2,j), [c(t1+j:end); zeros(j,1)], fc(:,3:5,j)]; 
    fn(:,:,j) = [fn(:,1:2,j), [n(t1+j:end); zeros(j,1)], fn(:,3:5,j)]; 
end
ft(ft==0)=NaN;
fy(fy==0)=NaN;
fc(fc==0)=NaN;
fn(fn==0)=NaN;
%% Forecast accuracy
load([pwd '\input\var.mat'])
nber(:, :, 1) = xlsread([pwd '\input\nberdates'],'nber1','E2:AZ1000');   
nber(:, :, 2) = xlsread([pwd '\input\nberdates'],'nber2','E2:AZ1000');  
nber(:, :, 3) = xlsread([pwd '\input\nberdates'],'nber3','E2:AZ1000');  
nber(:, :, 4) = xlsread([pwd '\input\nberdates'],'nber4','E2:AZ1000');
nber(:, :, 5) = xlsread([pwd '\input\nberdates'],'nber0','E2:AZ1000');
nber(nber==0)=NaN;
% Root mean square error (RMSE)
fcn7(filename, npr, fy, fc, fn, fy_var, fc_var, fn_var, nber,...
    dt, ft, ft_var)
% Mean absolute error (MAE)
fcn8(filename, npr, fy, fc, fn, fy_var, fc_var, fn_var, nber,...
    dt, ft, ft_var)
%Continous rank probability score (CRPS)
[crps_ecm1, crps_ecm2, crps_ecm3, crps_var1, crps_var2, crps_var3] =... 
    fcn9(filename, npr, nt, ny, L, fy, fc, fn, crps_mat, crps_var, nber);
%% Balanced growth ratios
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
    plot(xax1, zt(:, 1), 'k:', 'LineWidth', 1);   
    xlabel('Year')
    set(gca,'yLim',[0.1, 0.4],'YTick', 0.1 : 0.1 : 0.4)      
yyaxis right
    plot(xax1, zt(:, 2), 'b', 'LineWidth', 1);   
    set(gca,'yLim',[0.6, 0.9],'YTick', 0.6 : 0.1 : 0.9)  
    xlabel('Year')
    datetick('x','YYYY');
    set(gca, 'FontSize', 12)
    set(gca, 'FontName', 'Times New Roman')
    hBand = recessionplot('recessions', myrecession);   
    set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
    leg = legend({'$$\log (y_t / c_t)$$',...
        '$$\log [y_t/c_t \times (1-n_t) /n_t] $$'}, 'interpreter', 'latex');
    set(leg, 'Location', 'northwest')
    legend boxoff
saveas(gcf, [pwd '\output\balanced_growth_ratios.png'])
%% Laws of motion
figure
    subplot(3, 1, 1)
        plot(xax1, zt(:, 3), 'b', 'LineWidth', 1);  
        hline(0, 'k:')        
        xlabel('Year')
        set(gca,'yLim',[-0.03, 0.05],'YTick', -0.03 : 0.02 : 0.05)      
        datetick('x','YYYY');
        set(gca, 'FontSize', 12)
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot;   
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)       
        xlabel('Year')
        title('$$\Delta \log (y_t / n_t)$$', 'interpreter', 'latex')
    subplot(3, 1, 2)
        plot(xax1, zt(:, 4), 'b', 'LineWidth', 1); 
        hline(0, 'k:')        
        xlabel('Year')
        set(gca,'yLim',[-0.03, 0.05],'YTick', -0.03 : 0.02 : 0.05)      
        datetick('x','YYYY');
        set(gca, 'FontSize', 12)
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot;   
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)        
        xlabel('Year')
        title('$$\Delta \log (c_t / n_t)$$', 'interpreter', 'latex')        
    subplot(3, 1, 3)
        plot(xax1, zt(:, 5), 'b', 'LineWidth', 1);
        hline(0, 'k:')
        xlabel('Year')
        set(gca,'yLim', [-0.05, 0.07],'YTick', -0.05 : 0.03 : 0.07)      
        datetick('x','YYYY');
        set(gca, 'FontSize', 12)
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot;   
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)       
        xlabel('Year')
        title('$$\Delta \log (1-n_t / n_t)$$', 'interpreter', 'latex')        
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf, [pwd '\output\laws_of_motion'],'-dpdf')   
%% y, c, and n series
figure
	subplot(3, 1, 1)
        plot(xax, y / 10000, 'b', 'LineWidth', 1.5);
        axis([-Inf Inf 0.5 3.5])
        datetick('x','YYYY');
        set(gca, 'FontSize', 12) 
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot; 
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
        set(gca,'yLim',[0.5, 3.5],'YTick', 0.5 : 1 : 3.5)     
        xlabel('Year')
        ylabel('10,000 of chained 2012 dollars')  
        title('\it Real output per capita', 'FontWeight','Normal')
    subplot(3, 1, 2)	
        plot(xax, c / 10000, 'b', 'LineWidth', 1.5);
        axis([-Inf Inf 0.5 3.5])    
        datetick('x','YYYY');
        set(gca, 'FontSize', 12) 
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot;
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
        set(gca,'yLim',[0.5, 3.5],'YTick', 0.5 : 1 : 3.5)      
        xlabel('Year')
        ylabel('10,000 of chained 2012 dollars')
        title('\it Real consmuption per capita', 'FontWeight','Normal')     
    subplot(3, 1, 3)
        plot(xax, n * 100, 'b', 'LineWidth', 1.5);
        axis([-Inf Inf 32 41])
        datetick('x','YYYY');
        set(gca, 'FontSize', 12) 
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot;
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
        xlabel('Year') 
        ylabel('Percentage')
        title('\it Fraction of time spent working', 'FontWeight','Normal')        
        set(gca,'yLim',[32, 41],'YTick', 32 : 3 : 41)  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf, [pwd '\output\data_series'],'-dpdf')        
%% State variable series
figure
    plot(xax1, st(:, 1), 'k:', 'LineWidth', 0.5);
    hold on
    plot(xax1(3 : end), s(3 : end, 1), 'b', 'LineWidth', 1.5);
    hline(0, 'k:')
    axis([-Inf Inf -0.035 0.045])
    ax = gca;
    ax.XTick = xax;
    datetick('x','YYYY');
    set(gca, 'FontName', 'Times New Roman')
    hBand = recessionplot;
    set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)  
    set(gca, 'FontSize', 14) 
    xlabel('Year')
    set(gca, 'yLim', [-0.02, 0.03],'YTick', -0.02 : 0.01 : 0.03)       
    leg = legend({'Initial values', 'Fitted values'},...
        'interpreter', 'latex');
    set(leg, 'Location', 'southwest', 'FontSize', 12)
    legend boxoff
    saveas(gcf, [pwd '\output\state_variable_g.png'])   
figure
    plot(xax1, 1 ./ (exp(st(:, 2)) + 1), 'k:', 'LineWidth', 0.5);
    hold on
    plot(xax1(3 : end), 1 ./ (exp(s(3 : end, 2)) + 1) ,...
        'b', 'LineWidth', 1.5); 
    axis([-Inf Inf 0.37 0.47])    
    ax = gca;
    ax.XTick = xax;
    datetick('x','YYYY');
    set(gca, 'FontName', 'Times New Roman')
    hBand = recessionplot('recessions', myrecession);  
    set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)        
    xlabel('Year')     
    set(gca, 'FontSize', 14) 
    xlabel('Year')
    set(gca, 'yLim', [0.37, 0.47],'YTick', 0.37 : 0.02 : 0.47)  
    leg = legend({'Initial values', 'Fitted values'},...
        'interpreter', 'latex');
    set(leg, 'Location', 'southwest', 'FontSize', 12)
    legend boxoff    
    saveas(gcf, [pwd '\output\state_variable_d.png'])        
figure
    plot(xax1, st(:, 3), 'k:', 'LineWidth', 0.5);
    hold on
    plot(xax1(3 : end), s(3 : end, 3), 'b', 'LineWidth', 1.5);  
    axis([-Inf Inf 0.26 0.46])    
    ax = gca;
    ax.XTick = xax;
    datetick('x','YYYY');
    set(gca, 'yLim', [0.26, 0.46],'YTick', 0.26 : 0.04 : 0.46)      
    set(gca, 'FontSize', 14) 
    set(gca, 'FontName', 'Times New Roman')
    hBand = recessionplot('recessions', myrecession); 
    set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1) 
    xlabel('Year')      
    leg = legend({'Initial values', 'Fitted values'},...
        'interpreter', 'latex');
    set(leg, 'Location', 'northwest', 'FontSize', 12)
    legend boxoff    
    saveas(gcf, [pwd '\output\state_variable_a.png'])   
%% VAR parameter stability
b = cell2mat(out(:, 4));
ax1 = 1.96;
figure
k = 0;  
for ii = 1 : 7
    for j = 1 : 3   
        k = k + 1;
        subplot(7, 3, k)
        bj = b(ii : 7 : end, j);
        tj = b(ii : 7 : end, j + 3);
        sj = bj ./ tj;
        x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
        inBetween = [(bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
        fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
        hold on
        plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)
        hold on 
        plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
        hold on 
        plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5) 
        hline(0, 'k:')
        ax = gca;
        ax.XTick = xax;
        if ii == 1
            set(gca,'yLim',[-0.3, 0.3],'YTick', -0.3 : 0.3 : 0.3)
        elseif ii == 2
            set(gca,'yLim',[0, 10],'YTick', 0 : 5 : 10)          
        elseif ii == 3
            set(gca,'yLim',[-0.5, 1.5],'YTick', -0.5 : 1 : 1.5)    
        elseif ii == 4
            set(gca,'yLim',[-4, 0],'YTick', -4 : 2 : 0)     
        elseif ii == 5
            set(gca,'yLim',[-10, 0],'YTick', -10 : 5 : 0)    
        elseif ii == 6
            set(gca,'yLim',[-0.6, 0],'YTick', -0.6 : 0.3 : 0)  
        elseif ii == 7
            set(gca,'yLim',[0, 5],'YTick', 0 : 2.5 : 5)                
        end     
        if ii == 1 && j == 1
            title('[$\hat{s}_{t,1}\sim$const]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 1 && j == 2
            title('[$\hat{s}_{t,2}\sim$const]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 1 && j == 3
            title('[$\hat{s}_{t,3}\sim$const]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)    
        elseif ii == 2 && j == 1
            title('[$\hat{s}_{t,1}\sim\hat{s}_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 2 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 2 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 1
            title('[$\hat{s}_{t,1}\sim \hat{s}_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)      
        elseif ii == 4 && j == 1
            title('[$\hat{s}_{t,1}\sim \hat{s}_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 4 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 4 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)        
        elseif ii == 5 && j == 1
            title('[$\hat{s}_{t,1}\sim \hat{s}_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 5 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 5 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)             
        elseif ii == 6 && j == 1
            title('[$\hat{s}_{t,1}\sim \hat{s}_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 6 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 6 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)   
        elseif ii == 7 && j == 1
            title('[$\hat{s}_{t,1}\sim \hat{s}_{t-2,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 7 && j == 2
            title('[$\hat{s}_{t,2}\sim\hat{s}_{t-2,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 7 && j == 3
            title('[$\hat{s}_{t,3}\sim\hat{s}_{t-2,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)                 
        end        
        datetick('x','YYYY');  
        ax.XTickLabel = {'1980', '', '2000', '', '2020'};         
        set(gca, 'FontSize', 10)          
        hBand = recessionplot('recessions', myrecession); 
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)
        set(gca, 'FontName', 'Times New Roman')
    end       
end    
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf, [pwd '\output\var_process'],'-dpdf')   
%% ECM parameter stability
b = cell2mat(out(:, 11));
ax1 = 1.96;
figure
k = 0;   
for ii = 1 : 7
    for j = 1 : 3  
        k = k + 1;
        subplot(7, 3, k)
        bj = b(ii : 7 : end, j);
        tj = b(ii : 7 : end, j + 3);
        sj = bj ./ tj;
        x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
        inBetween = [(bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
        fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
        hold on
        plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)
        hold on 
        plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
        hold on 
        plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5)    
        hline(0, 'k:')
        ax = gca;
        ax.XTick = xax;
        if ii == 1
            set(gca,'yLim',[-0.003, 0.003],'YTick', -0.003 : 0.003 : 0.003)
        elseif ii == 2
            set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.15 : 0.15)          
        elseif ii == 3
            set(gca,'yLim',[-1, 2],'YTick', -1 : 1.5 : 2)    
        elseif ii == 4
            set(gca,'yLim',[-1, 1.5],'YTick', -1 : 1.25 : 1.5)     
        elseif ii == 5
            set(gca,'yLim',[-0.2, 0.8],'YTick', -0.2 : 0.5 : 0.8)    
        elseif ii == 6
            set(gca,'yLim',[-1, 6],'YTick', -1 : 3.5 : 6)  
        elseif ii == 7
            set(gca,'yLim',[-2, 2],'YTick', -2 : 2 : 2)                
        end   
        if ii == 1 && j == 1
            title('[$\Delta x_{t,1}^o\sim$const]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)
        elseif ii == 1 && j == 2
             title('[$\Delta x_{t,2}\sim$const]',...
             'FontName', 'Times New Roman', 'FontWeight', 'normal',...
             'Interpreter', 'latex')
        elseif ii == 1 && j == 3
            title('[$\Delta x_{t,3}^o\sim$const]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex') 
        elseif ii == 2 && j == 1
            title('[$\Delta x_{t,1}^o\sim r_{t-1,1}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 2 && j == 2
             title('[$\Delta x_{t,2}^o\sim r_{t-1,1}^o$]',...
             'FontName', 'Times New Roman',...
             'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 2 && j == 3
            title('[$\Delta x_{t,3}^o\sim r_{t-1,1}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')     
        elseif ii == 3 && j == 1
            title('[$\Delta x_{t,1}^o\sim \Delta x_{t-1,1}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 3 && j == 2
             title('[$\Delta x_{t,2}^o\sim \Delta x_{t-1,1}^o$]',...
             'FontName', 'Times New Roman',...
             'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 3 && j == 3
            title('[$\Delta x_{t,3}^o\sim \Delta x_{t-1,1}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')             
        elseif ii == 4 && j == 1
            title('[$\Delta x_{t,1}^o\sim \Delta x_{t-1,2}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 4 && j == 2
            title('[$\Delta x_{t,2}^o\sim \Delta x_{t-2,1}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 4 && j == 3
            title('[$\Delta x_{t,3}^o\sim \Delta x_{t-1,2}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')   
        elseif ii == 5 && j == 1
            title('[$\Delta x_{t,1}^o\sim \Delta x_{t-1,3}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 5 && j == 2
            title('[$\Delta x_{t,2}^o\sim \Delta x_{t-1,3}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 5 && j == 3
            title('[$\Delta x_{t,3}^o\sim \Delta x_{t-1,3}^o$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex') 
        elseif ii == 6 && j == 1
            title('[$\Delta x_{t,1}^o\sim \Delta \hat{s}_{t,1}$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 6 && j == 2
           title('[$\Delta x_{t,2}^o\sim \Delta \hat{s}_{t,1}$]',...
           'FontName', 'Times New Roman',...
           'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 6 && j == 3
           title('[$\Delta x_{t,3}^o\sim \Delta \hat{s}_{t,1}$]',...
           'FontName', 'Times New Roman',...
           'FontWeight', 'normal', 'Interpreter', 'latex')   
        elseif ii == 7 && j == 1
            title('[$\Delta x_{t,1}^o\sim \Delta \hat{s}_{t,3}$]',...
            'FontName', 'Times New Roman',...
            'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 7 && j == 2
           title('[$\Delta x_{t,2}^o\sim \Delta \hat{s}_{t,3}$]',...
           'FontName', 'Times New Roman',...
           'FontWeight', 'normal', 'Interpreter', 'latex')
        elseif ii == 7 && j == 3
           title('[$\Delta x_{t,3}^o\sim \Delta \hat{s}_{t,3}$]',...
           'FontName', 'Times New Roman',...
           'FontWeight', 'normal', 'Interpreter', 'latex')               
        end        
        datetick('x','YYYY');  
        ax.XTickLabel = {'1980', '', '2000', '', '2020'};         
        set(gca, 'FontSize', 10)          
        hBand = recessionplot('recessions', myrecession);
        set(hBand,'FaceColor', 'r', 'FaceAlpha', 0.1)
        set(gca, 'FontName', 'Times New Roman')
    end      
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf, [pwd '\output\ecm_process'],'-dpdf')   
%% ECM components
b = cell2mat(out(:, 11));
fig = figure;
for j = 1 : 3     
    subplot(3, 1, j)
        bj = b(2 : 7 : end, j);
        tj = b(2 : 7 : end, j + 3);
        sj = bj ./ tj;
        x2 = [xax(t1 : t2)',  fliplr(xax(t1 : t2)')];
        inBetween = [ (bj - ax1 * sj)', fliplr((bj + ax1 * sj)')];
        fill(x2, inBetween, [0.8 0.8157 1], 'LineStyle','none');
        hold on
        plot(xax(t1 : t2), bj, 'b', 'LineWidth', 1.5)  
        hold on
        plot(xax(t1 : t2), bj - ax1 * sj , '--b', 'LineWidth', 0.5)      
        hold on 
        plot(xax(t1 : t2), bj + ax1 * sj , '--b', 'LineWidth', 0.5)            
        hline(0, 'k:')
        ax = gca;
        ax.XTick = xax;
        yt = yticks;      
        datetick('x','YYYY');  
        set(gca, 'FontSize', 14)
        if j == 1
            set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)  
        elseif j == 2
            set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)   
        elseif j == 3
            set(gca,'yLim',[-0.15, 0.15],'YTick', -0.15 : 0.05 : 0.15)              
        end
        set(gca, 'FontName', 'Times New Roman')        
        hBand = recessionplot('recessions', myrecession); 
        set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)  
        xlabel('Year')  
        if j == 1
        title( '{\it Equilibrium correction coefficient for $$\Delta \ln(y_t/n_t)$$}',...
         'interpreter', 'latex')
        elseif j == 2
        title('{\it Equilibrium correction coefficient for $$\Delta \ln(c_t/n_t)$$}',...
         'interpreter', 'latex')
        elseif j == 3
        title('{\it Equilibrium correction coefficient for $$\Delta \ln((1-n_t)/n_t)$$}',...
         'interpreter', 'latex')
        end
end  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperType', 'A4');
size = get(gcf, 'PaperSize');
width = size(1);
height = size(2);
set(gcf, 'PaperPosition', [0, 0, width, height])
print(gcf, [pwd '\output\ecm_correction_coefficients'],'-dpdf')
%% Figures 1-4
for i = 1 : 3
    if i == 1
        z = fy / 10000;
        zvar = fy_var / 10000;
        v = 'Y';
        ftt = ft(:, 1) / 10000;
        dtt = dt(:, 1) / 10000;
        ztt = ft_var(:, 1) / 10000;        
    elseif i == 2
        z = fc / 10000;
        zvar = fc_var / 10000;        
        v = 'C';  
        ftt = ft(:, 2) / 10000;
        dtt = dt(:, 2) / 10000;
        ztt = ft_var(:, 2) / 10000;        
    else
        z = fn * 100;
        zvar = fn_var * 100;        
        v = 'N';          
        ftt = ft(:, 3) * 100;
        dtt = dt(:, 3) * 100;
        ztt = ft_var(:, 3) * 100;        
    end
    sD1 = datenum('01-01-1985');
    eD1 = datenum('10-01-2020');
    xaxl = linspace(sD1, eD1, 144)';

    for j = 1 : npr - 1
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
        h1 = plot(fxax(:, j), z(:, 3, j));
        hold on
        h2 = plot(fxax(:, j), z(:, 4, j));  
        hold on
        h4 = plot(fxax(:, j), z(:, 5, j), 'Color', 'k', 'LineStyle',...
            '--', 'LineWidth', 0.1);  
        hold on
        plot(fxax(:, j), z(:, 6, j),'Color', 'k', 'LineStyle', '--',...
                'LineWidth', 0.1); 
        hold on
        h3 = plot(fxax(:, j), zvar(:, 4, j));  
        hold on
        h5 = plot(fxax(:, j), zvar(:, 5, j), 'Color', [0.5 0.5 0.5],...
            'LineStyle', '-.', 'LineWidth', 0.1);   
        hold on
        plot(fxax(:, j), zvar(:, 6, j), 'Color', [0.5 0.5 0.5],...
            'LineStyle', '-.', 'LineWidth', 0.1);         
        hold on
        plot(xaxl, -100, 'w')
        set(h1, 'LineWidth', 2.5, 'Color', 'b')
        set(h2, 'LineWidth', 1, 'Color', 'k')
        set(h3, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])        
        if i == 1
            set(gca,'yLim',[1.4 3.2],'YTick', 1.4 : 0.36 : 3.2)  
            ylabel('10,000 of chained 2012 dollars')
        elseif i == 2 
            set(gca,'yLim',[1. 2.6],'YTick', 1 : 0.32 : 2.6)
            ylabel('10,000 of chained 2012 dollars')            
        elseif i == 3
            set(gca,'yLim',[30 40],'YTick', 30 : 2 : 40) 
            ylabel('Percentage')            
        end               
        ax = gca;
        ax.XTick = xax;
        ax.TitleFontSizeMultiplier = 1.5;        
        datetick('x','YYYY')
        xlim([datenum('01-01-1985') datenum('10-01-2020')])
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot('recessions', myrecession);
        set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
        xlabel('Year')     
        leg = legend([h1, h2, h4, h3, h5],...
            {'Data', 'VAR-ECM point forecast',...
            '95% confidence interval for VAR-ECM forecast',...
            'Baseline VAR point forecast',...
            '95% confidence interval for baseline VAR forecast'});
        legend boxoff
        set(leg, 'Location', 'Northwest')
        set(gca, 'FontName', 'Times New Roman')  
        set(gca, 'FontSize', 16)        
        saveas(gcf, [pwd '\output\f', num2str(i) num2str(j) '.png'])
    end
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
        h1 = plot(fxax0, dtt);
        hold on
        h2 = plot(fxax0, ftt);  
        hold on
        h3 = plot(fxax0, ztt);          
        hold on
        plot(xaxl, -100, 'w')
        set(h1, 'LineWidth', 2.5, 'Color', 'b')
        set(h2, 'LineWidth', 1, 'Color', 'k')
        set(h3, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])        
        if i == 1
            set(gca,'yLim',[1.4 3.2],'YTick', 1.4 : 0.36 : 3.2)  
            ylabel('10,000 of chained 2012 dollars')
        elseif i == 2 
            set(gca,'yLim',[1. 2.6],'YTick', 1 : 0.32 : 2.6)
            ylabel('10,000 of chained 2012 dollars')            
        elseif i == 3
            set(gca,'yLim',[30 40],'YTick', 30 : 2 : 40) 
            ylabel('Percentage')            
        end               
        ax = gca;
        ax.XTick = xax;
        ax.TitleFontSizeMultiplier = 1.5;        
        datetick('x','YYYY')
        xlim([datenum('01-01-1985') datenum('10-01-2020')])
        set(gca, 'FontName', 'Times New Roman')
        hBand = recessionplot('recessions', myrecession);
        set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
        xlabel('Year')     
        leg = legend([h1, h2, h3],...
            {'Data', 'VAR-ECM fitted value',...
            'Baseline VAR fitted value'});
        legend boxoff
        set(leg, 'Location', 'Northwest')
        set(gca, 'FontName', 'Times New Roman')  
        set(gca, 'FontSize', 16)        
        saveas(gcf, [pwd '\output\f', num2str(i) '0' '.png'])      
end
%% Figures recessions' windows
rec1(:, :, 1) = xlsread([pwd '\input\recession1'],'nber1','E2:Z1000');   
rec1(:, :, 2) = xlsread([pwd '\input\recession1'],'nber2','E2:Z1000');  
rec1(:, :, 3) = xlsread([pwd '\input\recession1'],'nber3','E2:Z1000');  
rec1(:, :, 4) = xlsread([pwd '\input\recession1'],'nber4','E2:Z1000');
rec1(:, :, 5) = xlsread([pwd '\input\recession1'],'nber0','E2:Z1000');  

rec2(:, :, 1) = xlsread([pwd '\input\recession2'],'nber1','E2:Z1000');   
rec2(:, :, 2) = xlsread([pwd '\input\recession2'],'nber2','E2:Z1000');  
rec2(:, :, 3) = xlsread([pwd '\input\recession2'],'nber3','E2:Z1000');  
rec2(:, :, 4) = xlsread([pwd '\input\recession2'],'nber4','E2:Z1000');
rec2(:, :, 5) = xlsread([pwd '\input\recession2'],'nber0','E2:Z1000');

rec3(:, :, 1) = xlsread([pwd '\input\recession3'],'nber1','E2:Z1000');   
rec3(:, :, 2) = xlsread([pwd '\input\recession3'],'nber2','E2:Z1000');  
rec3(:, :, 3) = xlsread([pwd '\input\recession3'],'nber3','E2:Z1000');  
rec3(:, :, 4) = xlsread([pwd '\input\recession3'],'nber4','E2:Z1000');
rec3(:, :, 5) = xlsread([pwd '\input\recession3'],'nber0','E2:Z1000');

for i = 1 : 3
    if i == 1
        z = fy / 10000;
        zvar = fy_var / 10000;
        v = 'Y';
        ftt = ft(:, 1) / 10000;
        dtt = dt(:, 1) / 10000;
        ztt = ft_var(:, 1) / 10000;
    elseif i == 2
        z = fc / 10000;
        zvar = fc_var / 10000;        
        v = 'C';        
        ftt = ft(:, 2) / 10000;
        dtt = dt(:, 2) / 10000;  
        ztt = ft_var(:, 2) / 10000;        
    else
        z = fn * 100;
        zvar = fn_var * 100;        
        v = 'N';      
        ftt = ft(:, 3) * 100;
        dtt = dt(:, 3) * 100;  
        ztt = ft_var(:, 3) * 100;        
    end    
    sD1 = datenum('01-01-1985');
    eD1 = datenum('10-01-2020');
    xaxl = linspace(sD1, eD1, 144)';
    for r = 1 : 3
        fig = figure('units','normalized','outerposition',[0 0 1 1]);       
        if r == 1
            rec = rec1;
        elseif r == 2 
            rec = rec2;
        elseif r == 3
            rec = rec3;
        end
        rec(rec==0)=NaN;
        j = 0;
        for s = 1 : 4
            subplot(2, 2, s)
            if s == 1
                h2 = plot(fxax0, ftt .* rec(:, 1, 5)); 
                hold on
                h3 = plot(fxax0, ztt .* rec(:, 1, 5)); 
                hold on
                h1 = plot(fxax0, dtt .* rec(:, 1, 5)  , '-');
                set(h1, 'LineWidth', 2.5, 'Color', 'b')
                set(h2, 'LineWidth', 1, 'Color', 'k')
                set(h3, 'LineWidth', 1, 'Color', [0.7 0.7 0.7])                
            else      
                j = j + 1;
                h2 = plot(fxax(:, j), z(:, 4, j) .* rec(:, 1, j));  
                hold on
                h4 = plot(fxax(:, j), z(:, 5, j) .* rec(:, 1, j),...
                    '--', 'Color', 'k', 'LineWidth', 0.5);  
                hold on
                plot(fxax(:, j), z(:, 6, j) .* rec(:, 1, j),...
                    '--', 'Color', 'k', 'LineWidth', 0.5);  
                hold on
                h3 = plot(fxax(:, j), zvar(:, 4, j) .* rec(:, 1, j));  
                hold on
                h5 = plot(fxax(:, j), zvar(:, 5, j) .* rec(:, 1, j),...
                    'Color', [0.5 0.5 0.5], 'LineStyle', '-.',...
                    'LineWidth', 0.1);   
                hold on
                plot(fxax(:, j), zvar(:, 6, j) .* rec(:, 1, j),...
                    'Color', [0.5 0.5 0.5], 'LineStyle', '-.',...
                    'LineWidth', 0.1);    
                hold on
                h1 = plot(fxax(:, j), z(:, 3, j) .* rec(:, 1, j)  , '-');
                set(h1, 'LineWidth', 2.5, 'Color', 'b')
                set(h2, 'LineWidth', 1, 'Color', 'k')
                set(h3, 'LineWidth', 1, 'Color', [0.5 0.5 0.5]) 
            end
            hold on
            plot(xaxl, -100, 'w')
            hold on
            if i == 1 && r == 1
                set(gca,'yLim', [1.5, 1.8])  
                ylabel('10,000 of chained 2012 dollars')
            elseif i == 1 && r == 2
                set(gca,'yLim', [1.9, 2.35])  
                ylabel('10,000 of chained 2012 dollars')
            elseif i == 1 && r == 3
                set(gca,'yLim', [2.1, 2.6])  
                ylabel('10,000 of chained 2012 dollars')                  
            elseif i == 2 && r == 1               
                set(gca,'yLim', [1.3, 1.5])
                ylabel('10,000 of chained 2012 dollars')      
            elseif i == 2 && r == 2               
                set(gca,'yLim', [1.55, 1.85])
                ylabel('10,000 of chained 2012 dollars')
            elseif i == 2 && r == 3               
                set(gca,'yLim', [1.8, 2.05])
                ylabel('10,000 of chained 2012 dollars')                   
            elseif i == 3 && r == 1               
                set(gca,'yLim', [32, 38])
                ylabel('Percentage') 
            elseif i == 3 && r == 2               
                set(gca,'yLim', [33, 39])
                ylabel('Percentage')            
            elseif i == 3 && r == 3               
                set(gca,'yLim', [31, 37])
                ylabel('Percentage')                   
            end                     
            hBand = recessionplot('recessions', myrecession);
            set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
            ax = gca;
            ax.XTick = xax;
            ax.TitleFontSizeMultiplier = 1.5;        
            datetick('x','QQ-YYYY')
            if r == 1 
                xlim([datenum('Q1-1990', 'QQ-YY')...
                      datenum('Q3-1992', 'QQ-YY')])
                xticks([datenum('Q1-1990', 'QQ-YY')...
                        datenum('Q3-1990', 'QQ-YY')...
                        datenum('Q1-1991', 'QQ-YY')...
                        datenum('Q3-1991', 'QQ-YY')...
                        datenum('Q1-1992', 'QQ-YY')...
                        datenum('Q3-1992', 'QQ-YY')])
                xticklabels({'Q1-1990', 'Q3-1990', 'Q1-1991', 'Q3-1991',...
                             'Q1-1992', 'Q3-1992'})                
            elseif r == 2 
                xlim([datenum('Q3-2000', 'QQ-YY')...
                      datenum('Q2-2003', 'QQ-YY')])
                xticks([datenum('Q3-2000', 'QQ-YY')...
                      datenum('Q1-2001', 'QQ-YY')...
                      datenum('Q3-2001', 'QQ-YY')...
                      datenum('Q1-2002', 'QQ-YY')...
                      datenum('Q3-2002', 'QQ-YY')...
                      datenum('Q1-2003', 'QQ-YY')])
                xticklabels({'Q3-2000', 'Q1-2001', 'Q3-2001', 'Q1-2002',...
                             'Q3-2002', 'Q1-2003'})                        
            elseif r == 3 
                xlim([datenum('Q2-2007', 'QQ-YY')...
                      datenum('Q4-2010', 'QQ-YY')])                
                xticks([datenum('Q2-2007', 'QQ-YY')...
                      datenum('Q4-2007', 'QQ-YY')...
                      datenum('Q2-2008', 'QQ-YY')...
                      datenum('Q4-2008', 'QQ-YY')...
                      datenum('Q2-2009', 'QQ-YY')...
                      datenum('Q4-2009', 'QQ-YY')...
                      datenum('Q2-2010', 'QQ-YY')...
                      datenum('Q4-2010', 'QQ-YY')])
                xticklabels({'Q2-2007', 'Q4-2007', 'Q2-2008', 'Q4-2008',...
                         'Q2-2009', 'Q4-2009', 'Q2-2010', 'Q4-2010'})                 
            end
            xlabel('Year')     
            set(gca, 'FontName', 'Times New Roman')  
            set(gca, 'FontSize', 10)   
            if s == 1
                leg = legend([h1, h2, h3],...
                     {'Data', 'VAR-ECM fitted value',...
                     'Baseline VAR fitted value'});           
                title('\it 0-step ahead', 'FontWeight','Normal')                 
            else
                leg = legend([h1, h2, h4, h3, h5],...
                     {'Data', 'VAR-ECM point forecast',...
                     '95% confidence interval for VAR-ECM forecast',...
                     'Baseline VAR point forecast',...
                     '95% confidence interval for baseline VAR forecast'});
                title(['\it' num2str(j) '-step ahead'],...
                    'FontWeight', 'Normal')                 
            end
            legend boxoff
            if i == 1 || i == 3
                set(leg, 'Location', 'Southwest', 'FontSize', 6)
            else
                set(leg, 'Location', 'Northwest', 'FontSize', 6)
            end
            saveas(gcf, [pwd '\output\rec_', num2str(r) 'var_'...
                num2str(i) '.png'])
        end
    end
end
%% Hedgehog graphs
figure
is = 0;    
for i = 1 : 3
     for r = 1 : 2
        is = is + 1;
        subplot(3, 2, is)
        for q = 1 : 2
            if q == 1              
                if i == 1
                    v1 = fy / 10000;
                    v2 = ft(:, 1) / 10000;
                    v3 = dt(:, 1) / 10000;
                elseif i == 2
                    v1 = fc / 10000;
                    v2 = ft(:, 2) / 10000;       
                    v3 = dt(:, 2) / 10000;
                elseif i == 3
                    v1 = fn * 100;
                    v2 = ft(:, 3) * 100;   
                    v3 = dt(:, 3) * 100;
                end    
            elseif q == 2
                if i == 1
                    v1 = fy_var / 10000;
                    v2 = ft_var(:, 1) / 10000;
                    v3 = dt(:, 1) / 10000;
                elseif i == 2
                    v1 = fc_var / 10000;
                    v2 = ft_var(:, 2) / 10000;       
                    v3 = dt(:, 2) / 10000;
                elseif i == 3
                    v1 = fn_var * 100;
                    v2 = ft_var(:, 3) * 100;   
                    v3 = dt(:, 3) * 100;
                end                
            end
            if r == 1
                dh = 20;
                J = 8;
            elseif r == 2
                dh = 62;
                J = 9;
            end        
            for j = 1 : J - 2
                jj = dh - 1 + j;              
                if q == 1
                    h2 = plot(fxax0(jj : jj + 3), [v2(jj, 1),...
                    v1(jj, 4, 1), v1(jj, 4, 2), v1(jj, 4, 3)]);
                    set(h2, 'Color', 'b', 'MarkerSize', 2.5,...
                        'Marker', 'o', 'LineWidth', 0.5)
                    hold on
                    scatter(fxax0(jj : jj), v2(jj, 1), 2.5, 'filled',...
                        'MarkerFaceColor', 'b')
                elseif q == 2
                    h3 = plot(fxax0(jj : jj + 3), [v2(jj, 1),...
                    v1(jj, 4, 1), v1(jj, 4, 2), v1(jj, 4, 3)]);
                    set(h3, 'Color', [0.7 0.7 0.7], 'MarkerSize', 2.5,...
                    'Marker', 'o', 'LineWidth', 0.5)
                    hold on
                    scatter(fxax0(jj : jj), v2(jj, 1), 2.5, 'filled',...
                        'MarkerFaceColor', [0.7 0.7 0.7])
                end
                hold on
            end
            h1 = plot(fxax0(dh : dh + J), v3(dh : dh + J)); 
            set(h1, 'Color', 'k', 'LineWidth', 1.5)
            xlabel('Year')                
            if i == 1 && r == 1
                set(gca,'yLim', [1.6, 1.75])  
                ylabel('10,000 of chained 2012 dollars',...
                    'FontName', 'Times New Roman')  
                title('{\it Recession 1990-91: y}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')  
            elseif i == 1 && r == 2
                set(gca,'yLim', [2.1, 2.25])  
                ylabel('10,000 of chained 2012 dollars',...
                    'FontName', 'Times New Roman')  
                title('{\it Recession 2001: y}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')                                 
            elseif i == 2 && r == 1               
                set(gca,'yLim', [1.35, 1.41])
                yticks([1.35, 1.37, 1.39, 1.41])
                ylabel('10,000 of chained 2012 dollars',...
                    'FontName', 'Times New Roman')  
                title('{\it Recession 1990-91: c}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')                  
            elseif i == 2 && r == 2               
                set(gca,'yLim', [1.65, 1.74])
                yticks([1.65, 1.68, 1.71, 1.74])
                ylabel('10,000 of chained  2012 dollars')
                title('{\it Recession 2001: c}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')                                   
            elseif i == 3 && r == 1               
                 set(gca,'yLim', [33.5, 36.8])
                 %yticks([34, 35, 36])
                 title('{\it Recession 1990-91: n}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')                   
                 ylabel('Percentage',...
                    'FontName', 'Times New Roman')   
            elseif i == 3 && r == 2               
                 set(gca,'yLim', [34.5, 37.5])
                 yticks([34.5, 35.5, 36.5, 37.5])
                 ylabel('Percentage',...
                    'FontName', 'Times New Roman')       
                 title('{\it Recession 2001: n}',...
                    'FontWeight','Normal',...
                    'FontName', 'Times New Roman')                                              
            end                             
            hBand = recessionplot('recessions', myrecession);
            set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
            ax = gca;
            ax.XTick = xax;  
            ax.TitleFontSizeMultiplier = 1.5; 
            set(gca, 'FontName', 'Times New Roman')
            datetick('x','QQ-YY', 'keepticks')
            if r == 1 
                    xlim([datenum('Q4-1989', 'QQ-YY')...
                          datenum('Q3-1992', 'QQ-YY')])
                    xticks([datenum('Q1-1990', 'QQ-YY')...
                            datenum('Q2-1992', 'QQ-YY')])
                    xticklabels({'Q1-1990', 'Q2-1992'})   
            elseif r == 2 
                    xlim([datenum('Q2-2000', 'QQ-YY')...
                          datenum('Q1-2003', 'QQ-YY')])
                    xticks([datenum('Q3-2000', 'QQ-YY')...
                          datenum('Q4-2002', 'QQ-YY')])
                    xticklabels({'Q3-2000', 'Q4-2002'})                   
            end       
            axis square;
        end
        leg = legend([h1, h2, h3], {'Data', 'VAR-ECM', 'Benchmark VAR'});
        legend boxoff
        set(leg, 'Location', 'Northwest', 'FontSize', 7)
        set(gca, 'FontName', 'Times New Roman')  
        set(gca, 'FontSize', 9)         
    end
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperType', 'A4');
    size = get(gcf, 'PaperSize');
    width = size(1);
    height = size(2);
    set(gcf, 'PaperPosition', [0, 0, width, height])
	print(gcf, [pwd '\output\hedgehog_1990-2001'],'-dpdf')
end
% 2007-09 Great Recession
figure
is = 0; 
r = 3;
dh = 89;
J = 12; 
for i = 1 : 3
    is = is + 1;
    subplot(3, 2, is)
    for q = 1 : 2
        if q == 1              
            if i == 1
                v1 = fy / 10000;
                v2 = ft(:, 1) / 10000;
                v3 = dt(:, 1) / 10000;
            elseif i == 2
                v1 = fc / 10000;
                v2 = ft(:, 2) / 10000;       
                v3 = dt(:, 2) / 10000;
            elseif i == 3
                v1 = fn * 100;
                v2 = ft(:, 3) * 100;   
                v3 = dt(:, 3) * 100;
            end    
       elseif q == 2
            if i == 1
                v1 = fy_var / 10000;
                v2 = ft_var(:, 1) / 10000;
                v3 = dt(:, 1) / 10000;
            elseif i == 2
                v1 = fc_var / 10000;
                v2 = ft_var(:, 2) / 10000;       
                v3 = dt(:, 2) / 10000;
            elseif i == 3
                v1 = fn_var * 100;
                v2 = ft_var(:, 3) * 100;   
                v3 = dt(:, 3) * 100;
             end                
        end
        for j = 1 : J - 2
            jj = dh - 1 + j;   
            if q == 1
                scatter(fxax0(jj), v2(jj, 1), 2.5, 'filled',...
                'MarkerFaceColor', 'b')              
            elseif q == 2
               scatter(fxax0(jj), v2(jj, 1), 2.5, 'filled',...
               'MarkerFaceColor', [0.7 0.7 0.7])  
            end
            hold on
            if q == 1
                h2 = plot(fxax0(jj : jj + 3), [v2(jj, 1),...
                v1(jj, 4, 1), v1(jj, 4, 2), v1(jj, 4, 3)]);
                set(h2, 'Color', 'b', 'MarkerSize', 2.5,...
                'Marker', 'o', 'LineWidth', 0.5)          
            elseif q == 2
                h3 = plot(fxax0(jj : jj + 3), [v2(jj, 1),...
                v1(jj, 4, 1), v1(jj, 4, 2), v1(jj, 4, 3)]);
                set(h3, 'Color', [0.7 0.7 0.7], 'MarkerSize', 2.5,...
                'Marker', 'o', 'LineWidth', 0.5)
            end
            hold on
         end
         h1 = plot(fxax0(dh : dh + J), v3(dh : dh + J)); 
         set(h1, 'Color', 'k', 'LineWidth', 1.5)
         xlabel('Year')                                
         if i == 1 && r == 3
            set(gca,'yLim', [2.20, 2.53])  
            yticks([2.20, 2.31, 2.42, 2.53])                
            ylabel('10,000 of chained 2012 dollars',...
            'FontName', 'Times New Roman')    
            title('{\it Recession 2007-09: y}',...
            'FontWeight','Normal', 'FontName', 'Times New Roman')                                    
         elseif i == 2 && r == 3               
            set(gca,'yLim', [1.87, 1.96])
            yticks([1.87 1.90, 1.93, 1.96])
            ylabel('10,000 of chained 2012 dollars',...
            'FontName', 'Times New Roman')     
            title('{\it Recession 2007-09: c}',...
            'FontWeight','Normal','FontName', 'Times New Roman')                                       
         elseif i == 3 && r == 3               
            set(gca,'yLim', [32, 38])
            yticks([32, 34, 36, 38])
            ylabel('Percentage', 'FontName', 'Times New Roman')     
            title('{\it Recession 2007-09: n}',...
            'FontWeight','Normal', 'FontName', 'Times New Roman')                         
         end   
         hBand = recessionplot('recessions', myrecession);
         set(hBand, 'FaceColor', 'r', 'FaceAlpha', 0.1)
         ax = gca;
         ax.XTick = xax;  
         ax.TitleFontSizeMultiplier = 1.5; 
         set(gca, 'FontName', 'Times New Roman')
         datetick('x','QQ-YY', 'keepticks')
         xlim([datenum('Q1-2007', 'QQ-YY') datenum('Q4-2010', 'QQ-YY')])                
         xticks([datenum('Q2-2007', 'QQ-YY') datenum('Q3-2010', 'QQ-YY')])
         xticklabels({'Q2-2007', 'Q3-2010'})                 
         axis square;
    end
    leg = legend([h1, h2, h3], {'Data', 'VAR-ECM', 'Benchmark VAR'});
    legend boxoff
    set(leg, 'Location', 'Northwest', 'FontSize', 7)
    set(gca, 'FontName', 'Times New Roman')  
    set(gca, 'FontSize', 9)         
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperType', 'A4');
    size = get(gcf, 'PaperSize');
    width = size(1);
    height = size(2);
    set(gcf, 'PaperPosition', [0, 0, width, height])
	print(gcf, [pwd '\output\hedgehog_2007-2009'],'-dpdf')
end


