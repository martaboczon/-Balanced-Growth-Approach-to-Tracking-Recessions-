% FUNCTION: MEAN ABSOLUTE ERROR
function fcn8(filename, npr, fy, fc, fn, fy_var, fc_var, fn_var, nber,...
    dt, ft, ft_var)
K = 3; %Number of windows of analysis
j1 = 1;
for i = 1 : K
    j2 = j1 + 1;
    j3 = j1 + 2;
    j4 = j1 + 3;
    j5 = j1 + 4;
    xlswrite(filename, {['ECM: ' num2str(i)]},...
        'fstat-mae', ['A' num2str(j1)])      
    xlswrite(filename, {['VAR: ' num2str(i)]},...
        'fstat-mae', ['F' num2str(j1)])    
    xlswrite(filename, {'y'}, 'fstat-mae', ['B' num2str(j1)])
    xlswrite(filename, {'c'}, 'fstat-mae', ['C' num2str(j1)])
    xlswrite(filename, {'n'}, 'fstat-mae', ['D' num2str(j1)])
    xlswrite(filename, {'y'}, 'fstat-mae', ['G' num2str(j1)])
    xlswrite(filename, {'c'}, 'fstat-mae', ['H' num2str(j1)])
    xlswrite(filename, {'n'}, 'fstat-mae', ['I' num2str(j1)])
    xlswrite(filename, {'h=0'}, 'fstat-mae', ['A' num2str(j2)])
    xlswrite(filename, {'h=1'}, 'fstat-mae', ['A' num2str(j3)])
    xlswrite(filename, {'h=2'}, 'fstat-mae', ['A' num2str(j4)])
    xlswrite(filename, {'h=3'}, 'fstat-mae', ['A' num2str(j5)])
    xlswrite(filename, {'h=0'}, 'fstat-mae', ['F' num2str(j2)])
    xlswrite(filename, {'h=1'}, 'fstat-mae', ['F' num2str(j3)])
    xlswrite(filename, {'h=2'}, 'fstat-mae', ['F' num2str(j4)])
    xlswrite(filename, {'h=3'}, 'fstat-mae', ['F' num2str(j5)])    
    j1 = j1 + 6;
end
i = 1;
for w = 1 : K
    i = 1 + i;    
    xlswrite(filename, mean(abs((dt(:, 1) - ft(:, 1))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['B' num2str(i)]);   
    xlswrite(filename, mean(abs((dt(:, 2) - ft(:, 2))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['C' num2str(i)]);   
    xlswrite(filename, mean(abs((dt(:, 3) - ft(:, 3))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['D' num2str(i)]);    
    xlswrite(filename, mean(abs((dt(:, 1) - ft_var(:, 1))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['G' num2str(i)]);  
    xlswrite(filename, mean(abs((dt(:, 2) - ft_var(:, 2))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['H' num2str(i)]);    
    xlswrite(filename, mean(abs((dt(:, 3) - ft_var(:, 3))...
            .* nber(:, w, 5)), 'omitnan'), 'fstat-mae', ['I' num2str(i)]);         
    for j = 1 : npr - 1
        i = i + 1;
        % Mean Absolute Error (MAE)
        xlswrite(filename, mean(abs((fy(:, 3, j) - fy(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['B' num2str(i)]);   
        xlswrite(filename, mean(abs((fc(:, 3, j) - fc(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['C' num2str(i)]);   
        xlswrite(filename, mean(abs((fn(:, 3, j) - fn(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['D' num2str(i)]);    
        xlswrite(filename, mean(abs((fy_var(:, 3, j) - fy_var(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['G' num2str(i)]);  
        xlswrite(filename, mean(abs((fc_var(:, 3, j) - fc_var(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['H' num2str(i)]);    
        xlswrite(filename, mean(abs((fn_var(:, 3, j) - fn_var(:, 4, j))...
            .* nber(:, w, j)), 'omitnan'), 'fstat-mae', ['I' num2str(i)]);     
    end
    i = i + 2;
end