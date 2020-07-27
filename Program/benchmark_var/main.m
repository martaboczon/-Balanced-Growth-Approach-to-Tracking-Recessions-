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
%% Initial parameters
warning('off','MATLAB:xlswrite:AddSheet')
ts = 1.96;  % t-statistic threshold
npr = 4;    % foresacting horizon
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
startDate = datenum('01-01-1948');
endDate = datenum('04-01-2019');
xax = linspace(startDate, endDate, size(y, 1))';
 
startDate1 = datenum('04-01-1948');
endDate1 = datenum('04-01-2019');
xax1 = linspace(startDate1, endDate1, size(y, 1) - 1)';
 
sD1 = datenum('07-01-1985');
eD1 = datenum('07-01-2019');
fxax1 = linspace(sD1, eD1, nt)';
 
sD2 = datenum('10-01-1985');
eD2 = datenum('10-01-2019');
fxax2 = linspace(sD2, eD2, nt)';
 
sD3 = datenum('01-01-1986');
eD3 = datenum('01-01-2020');
fxax3 = linspace(sD3, eD3, nt)';
 
sD4 = datenum('04-01-1986');
eD4 = datenum('04-01-2020');
fxax4 = linspace(sD4, eD4, nt)';
fxax = [fxax1, fxax2, fxax3, fxax4];
%% Recursive estimation
for i = 1 : nt
%     pcav = pca([c(2 : t(i)) ./ c(1 : t(i) - 1),...
%         y(2 : t(i)) ./ y(1 : t(i) - 1)]);
%     c1 = pcav(1, 1);
%     c2 = pcav(2, 1);
%     au = c1 + c2;
%     c1 = c1 / au;
%     c2 = c2 / au;
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
    st = [zt3, zt4, zt5];
    cp = cp(2 : end, :);
    T0 = size(st, 1);
    % VAR process
    Y = st(5 : end, :) - st(4 : end - 1, :);
    X = [ones(T0 - 4, 1), st(4 : end - 1, :) - st(3 : end - 2, :),...
        st(3 : end - 2, :)]; 
    ny = size(Y, 2);
    nx = size(X, 2);
    [USURE_pi, USURE_sigma, USURE_std] = fcn1(X, Y, nx, ny);
    out{i, 4} = [USURE_pi, USURE_std];
    out{i, 5} = USURE_sigma;
    out{i, 6} = []; 
    % Tracking
    ft(i, :) =  fcn2(T0, npr, cp, zt, USURE_pi);
    % Forecasting
    [f, crps_mat{i}] =  fcn3(T0, ny, npr, L, cp, zt, USURE_pi, USURE_sigma); 
    for j = 1 : npr     
        fy(i, 1 : 5, j) = [i, T0 + j, f(j,:,1)]; 
        fc(i, 1 : 5, j) = [i, T0 + j, f(j,:,2)]; 
        fn(i, 1 : 5, j) = [i, T0 + j, f(j,:,3)];              
    end
end
for j = 1 : npr    
    fy(:,:,j) = [fy(:,1:2,j), [y(t1+j:end); zeros(j,1)], fy(:,3:5,j)]; 
    fc(:,:,j) = [fc(:,1:2,j), [c(t1+j:end); zeros(j,1)], fc(:,3:5,j)]; 
    fn(:,:,j) = [fn(:,1:2,j), [n(t1+j:end); zeros(j,1)], fn(:,3:5,j)]; 
end
ft(ft==0)=NaN;
fy(fy==0)=NaN;
fc(fc==0)=NaN;
fn(fn==0)=NaN;
for j = 1 : npr
    xlswrite(filename, fy(:, :, j), ['fy' num2str(j)], 'A2')
    xlswrite(filename, fc(:, :, j), ['fc' num2str(j)], 'A2')
    xlswrite(filename, fn(:, :, j), ['fn' num2str(j)], 'A2')      
end 
%% Parameter stability
% VAR
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
            set(gca,'yLim',[0, 0.01],'YTick', 0 : 0.005 : 0.01) 
        elseif ii == 2
            set(gca,'yLim',[-1.4, 0.2],'YTick', -1.4 : 0.8 : 0.2)   
        elseif ii == 3
            set(gca,'yLim',[-1.6, 0.8],'YTick', -1.6 : 1.2 : 0.8) 
        elseif ii == 4
            set(gca,'yLim',[-0.7, 0.7],'YTick', -0.7 : 0.7 : 0.7)      
        elseif ii == 5
            set(gca,'yLim',[-1.4, 0.2],'YTick', -1.4 : 0.8 : 0.2)
        elseif ii == 6
            set(gca,'yLim',[-1.4, 1],'YTick', -1.4 : 1.2 : 1)  
        elseif ii == 7
            set(gca,'yLim',[-0.6, 0.6],'YTick', -0.6 : 0.6 : 0)              
        end                    
        if ii == 1 && j == 1
            title('[$\Delta{x}_{t,1}\sim const$]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 1 && j == 2
            title('[$\Delta{x}_{t,2}\sim const$]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)
        elseif ii == 1 && j == 3
            title('[$\Delta{x}_{t,3}\sim const$]', 'FontName',...
            'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)      
        elseif ii == 2 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 2 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 2 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-1,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 3 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-1,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)       
        elseif ii == 4 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 4 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 4 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-1,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)    
        elseif ii == 5 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 5 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 5 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-2,1}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)   
        elseif ii == 6 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 6 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 6 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-2,2}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)  
        elseif ii == 7 && j == 1
            title('[$\Delta{x}_{t,1}\sim\Delta x_{t-2,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6) 
        elseif ii == 7 && j == 2
            title('[$\Delta{x}_{t,2}\sim\Delta x_{t-2,3}$]',...
             'FontName', 'Times New Roman', 'FontWeight', 'normal',...
             'Interpreter', 'latex', 'FontSize', 6) 
         elseif ii == 7 && j == 3
            title('[$\Delta{x}_{t,3}\sim\Delta x_{t-2,3}$]',...
            'FontName', 'Times New Roman', 'FontWeight', 'normal',...
            'Interpreter', 'latex', 'FontSize', 6)                      
        end
        datetick('x','YYYY');  
        ax.XTickLabel = {'1980', '', '2000', '', '2020'};         
        set(gca, 'FontSize', 10)          
        hBand = recessionplot;  
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
print(gcf, [pwd '\output\parameters'],'-dpdf')
ft_var = ft;
fc_var = fc;
fn_var = fn;
fy_var = fy;
crps_var = crps_mat;
save([pwd '\output\var.mat'],'fy_var', 'fc_var', 'fn_var', 'ft_var',...
'crps_var')
