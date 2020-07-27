% FUNCTION: CONTINOUS RANK PROBABILITY SCORE
function [crps_ecm1, crps_ecm2, crps_ecm3,...
    crps_var1, crps_var2, crps_var3] = ...
    fcn9(filename, npr, nt, ny, L, fy, fc, fn, crps_mat,...
    crps_var, nber)
K = 3; %Number of windows of analysis
j1 = 1;
for i = 1 : K
    j2 = j1 + 1;
    j3 = j1 + 2;
    j4 = j1 + 3;
    xlswrite(filename, {['ECM: ' num2str(i)]},...
        'fstat-crps', ['A' num2str(j1)])      
    xlswrite(filename, {['VAR: ' num2str(i)]},...
        'fstat-crps', ['F' num2str(j1)])    
    xlswrite(filename, {'y'}, 'fstat-crps', ['B' num2str(j1)])
    xlswrite(filename, {'c'}, 'fstat-crps', ['C' num2str(j1)])
    xlswrite(filename, {'n'}, 'fstat-crps', ['D' num2str(j1)])
    xlswrite(filename, {'y'}, 'fstat-crps', ['G' num2str(j1)])
    xlswrite(filename, {'c'}, 'fstat-crps', ['H' num2str(j1)])
    xlswrite(filename, {'n'}, 'fstat-crps', ['I' num2str(j1)])
    xlswrite(filename, {'h=1'}, 'fstat-crps', ['A' num2str(j2)])
    xlswrite(filename, {'h=2'}, 'fstat-crps', ['A' num2str(j3)])
    xlswrite(filename, {'h=3'}, 'fstat-crps', ['A' num2str(j4)])
    xlswrite(filename, {'h=1'}, 'fstat-crps', ['F' num2str(j2)])
    xlswrite(filename, {'h=2'}, 'fstat-crps', ['F' num2str(j3)])
    xlswrite(filename, {'h=3'}, 'fstat-crps', ['F' num2str(j4)])    
    j1 = j1 + 5;
end
crps_ecm1 = cell(2, 1);
crps_ecm2 = cell(2, 1);
crps_ecm3 = cell(2, 1);
crps_var1 = cell(2, 1);
crps_var2 = cell(2, 1);
crps_var3 = cell(2, 1);
for i = 1 : K
    crps_ecm1{i} = cell(npr - 1, ny);
    crps_ecm2{i} = cell(npr - 1, ny);
    crps_ecm3{i} = cell(npr - 1, ny);
    crps_var1{i} = cell(npr - 1, ny);
    crps_var2{i} = cell(npr - 1, ny);
    crps_var3{i} = cell(npr - 1, ny);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR-ECM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
step1 = cell2mat(crps_mat);
step2 = reshape(step1, nt, npr, L, ny);
step5 = zeros(npr - 1, ny);
for w = 1 : K       
    for j = 1 : npr - 1
        i = i + 1;
        for v = 1 : ny
            step3 = step2(:, j, :, v) ;
            step4 = reshape(step3, nt, L) .* nber(:, w, j);
            if v == 1
               step5(v, j) = crps(step4, fy(:, 3, j) .* nber(:, w, j));
               xlswrite(filename, step5(v, j), 'fstat-crps',...
                   ['B' num2str(i)]);   
               [crps_ecm1{w}{j, v}, crps_ecm2{w}{j, v},...
                   crps_ecm3{w}{j, v}] =...
                   crps(step4, fy(:, 3, j) .* nber(:, w, j));
            elseif v == 2               
               step5(v, j) = crps(step4, fc(:, 3, j) .* nber(:, w, j));
               xlswrite(filename, step5(v, j), 'fstat-crps',...
                   ['C' num2str(i)]); 
               [crps_ecm1{w}{j, v}, crps_ecm2{w}{j, v},...
                   crps_ecm3{w}{j, v}] =...
                   crps(step4, fc(:, 3, j) .* nber(:, w, j));               
            elseif v == 3               
              step5(v, j) = crps(step4, fn(:, 3, j) .* nber(:, w, j));
              xlswrite(filename, step5(v, j), 'fstat-crps',...
                  ['D' num2str(i)]);   
              [crps_ecm1{w}{j, v}, crps_ecm2{w}{j, v},...
                  crps_ecm3{w}{j, v}] =...
                  crps(step4, fn(:, 3, j) .* nber(:, w, j));              
            end
           
        end
    end
    i = i + 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
step1 = cell2mat(crps_var);
step2 = reshape(step1, nt, npr, L, ny);
step5 = zeros(npr - 1, ny);
for w = 1 : K       
    for j = 1 : npr - 1
        i = i + 1;
        for v = 1 : ny
            step3 = step2(:, j, :, v) ;
            step4 = reshape(step3, nt, L) .* nber(:, w, j);
            if v == 1
               step5(v, j) = crps(step4, fy(:, 3, j) .* nber(:, w, j));
               xlswrite(filename, step5(v, j), 'fstat-crps',...
                   ['G' num2str(i)]);   
               [crps_var1{w}{j, v}, crps_var2{w}{j, v},...
                   crps_var3{w}{j, v}] =...
                   crps(step4, fy(:, 3, j) .* nber(:, w, j));               
            elseif v == 2
               step5(v, j) = crps(step4, fc(:, 3, j) .* nber(:, w, j));
               xlswrite(filename, step5(v, j), 'fstat-crps',...
                   ['H' num2str(i)]); 
               [crps_var1{w}{j, v}, crps_var2{w}{j, v},...
                   crps_var3{w}{j, v}] =...
                   crps(step4, fc(:, 3, j) .* nber(:, w, j));                    
            elseif v == 3
               step5(v, j) = crps(step4, fn(:, 3, j) .* nber(:, w, j));
              xlswrite(filename, step5(v, j), 'fstat-crps',...
                  ['I' num2str(i)]);   
               [crps_var1{w}{j, v}, crps_var2{w}{j, v},...
                   crps_var3{w}{j, v}] =...
                   crps(step4, fn(:, 3, j) .* nber(:, w, j));                   
            end
        end
    end
    i = i + 2;
end