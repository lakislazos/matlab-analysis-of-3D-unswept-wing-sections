% house-keeping
clc; 
clear; 

% load airfoil data
file = input('Enter the airfoil data file: ', 's'); 
airfoil_data = load(file); % load the airfoil data from the specified file
fprintf('Airfoil data loaded: %d data points\n', size(airfoil_data,1)); % confirm data is loaded and print a message for the user

% input wing data
S = input('Enter the wing area in meters squared: '); 
taper_ratio = input('Enter the taper ratio: ');
AR = input('Enter the aspect ratio: ');
e = input('Enter the tip twist angle in degrees: ');

b = sqrt(AR * S); % calculate wing span 
c_root = 2*S / (b * (1 + taper_ratio)); % calculate root chord 
c_tip  = taper_ratio * c_root; % calculate tip chord

% input flight conditions
h = input('Enter the altitude in meters: ');
v_free = input('Enter the freestream velocity in meters/second: ');
[~, ~, ~, rho] = atmosisa(h); % air density at the specified altitude through atmosisa function

% extract columns from the imported file
AoA_data = airfoil_data(:, 1); % extract angle of attack in degrees from the file
AoA_data = deg2rad(AoA_data); % convert angles of attack to radians
Cl_data = airfoil_data(:, 2); % extract the sectional lift coefficient from the file
Cd_data = airfoil_data(:, 3); % extract the sectional drag coefficient from the file

% print a message for the user describing all available options
fprintf('\nSelect the Mode:\n')
fprintf('Mode 1: Fixed angle of attack and number of panels\n')
fprintf('Mode 2: Fixed number of panels and range of angles of attack from -20 to 30 degrees\n')
fprintf('Mode 3: Fixed angle of attack and range of number of panels (N = 10, 20, 50, 100)\n\n')

mode = input('Enter your desired mode as described above: ');

% set up AoA and N arrays based on the user's preference
if mode == 1
    AoA_infinity = input('Enter the angle of attack in degrees: ');
    N = input('Enter the number of panels: ');
    AoA_array = deg2rad(AoA_infinity); % convert the angle of attack to radians
    N_array = N; % store the N value in an array for consistency with the looping approach
elseif mode == 2
    N = input('Enter the number of panels: ');
    AoA_array = deg2rad(-20:1:30); % provide an array of angles of attack, converting them to radians
    N_array = N; % store the N value in an array for consistency with the looping approach
elseif mode == 3
    AoA_infinity = input('Enter the angle of attack in degrees: ');
    AoA_array = deg2rad(AoA_infinity); % convert the angle of attack to radians
    N_array = [10, 20, 50, 100]; % provide an array of numbers of panels
else
    error('Incorrect mode selected; enter 1, 2 or 3 for mode!') % produce an error for an incorrect user input
end

% preallocate matrices for results
CL_results = zeros(length(N_array), length(AoA_array)); % matrix to store C_L for each N and AoA entries
CD_results = zeros(length(N_array), length(AoA_array)); % matrix to store C_D for each N and AoA entries

% text file for final results
if mode == 1
    file_name = sprintf('results_for_mode=%d_AoA=%.1fdeg_N=%d.txt', mode, AoA_infinity, N);
elseif mode == 2
    file_name = sprintf('results_for_mode=%d_N=%d.txt', mode, N);
elseif mode == 3
    file_name = sprintf('results_for_mode=%d_AoA=%.1fdeg.txt', mode, AoA_infinity);
end

results_file = fopen(file_name, 'w'); % open the file using the fopen function

fprintf(results_file, '  3D Wing Analysis\n'); % create a title at the top
fprintf(results_file, '----------------------\n\n'); % create a separator for better readability

fprintf(results_file, '--- Flight Conditions ---\n'); % create a section title
fprintf(results_file, '  Altitude:            %.4f m\n', h);
fprintf(results_file, '  Freestream velocity: %.4f m/s\n', v_free);
fprintf(results_file, '  Air density:         %.4f kg/m^3\n', rho);

fprintf(results_file, '\n--- Wing Geometry ---\n'); 
fprintf(results_file, '  Wing area:           %.4f m^2\n', S);
fprintf(results_file, '  Taper ratio:         %.4f\n', taper_ratio);
fprintf(results_file, '  Aspect ratio:        %.4f\n', AR);
fprintf(results_file, '  Tip twist:           %.4f deg\n', e);
fprintf(results_file, '  Span:                %.4f m\n', b);
fprintf(results_file, '  Root chord:          %.4f m\n', c_root);
fprintf(results_file, '  Tip chord:           %.4f m\n', c_tip);

fprintf(results_file, '\n--- Flight Mode ---\n'); 
fprintf(results_file, '  Mode:                %d\n', mode);

fprintf(results_file, '\n--- Results ---\n');
fprintf(results_file, '  %-10s %-10s %-10s %-10s\n', 'AoA(deg)', 'N', 'CL', 'CD');
fprintf(results_file, '  %-10s %-10s %-10s %-10s\n', '--------', '--', '--', '--');

% set up convergence settings
tol = 1e-8; % convergence tolerance
max_iterations = 5000; % maximum number of iterations


for i = 1:length(N_array) % start the main loop (looping through the number of panels)
    N = N_array(i); % retrieve the current number of panels
    for j = 1:length(AoA_array) % start the secondary loop (looping through angles of attack)
        % discretize the wing
        delta_y = b / N; % width of 1 panel 
        y_i = -b/2 + delta_y/2 : delta_y : b/2 - delta_y/2; % center of the panels positions 
        c_i = c_root + (c_tip - c_root) * (2*abs(y_i)/b); % chord at centers of the panels
        e_y = deg2rad(e * (2*abs(y_i)/b)); % twist angle at centers of panels converted to radians

        G_i = (4*v_free*S*AoA_array(j)/b) * sqrt(1 - (y_i/(b/2)).^2); % initial elliptic circulation distribution
        
        % change the value of the damping factor depending on the number of panels; I decided to implement it this way after running the program multiple times and seeing that at higher numbers of panels, you need a lower damping factor to get convergence 
        if N <= 20 
            D = 0.07; 
        elseif N <= 50
            D = 0.04;
        else
            D = 0.01;
        end

        converged = false; % initialize a flag for convergence
        for iterations = 1:max_iterations % start the solver 

            % calculate the downwash at each panel using hshoe()
            W_i = zeros(1, N); % initialise downwash velocity array
            y1 = [y_i - delta_y/2; y_i + delta_y/2]'; % matrix of edge positions of each panel
            for k = 1:N
                W_i(k) = sum(hshoe(y_i(k), y1, G_i')); % sum all downwash contributions by using hshoe for each panel 
            end

            AoA_effective = e_y + atan((v_free*sin(AoA_array(j)) - W_i)./(v_free*cos(AoA_array(j)))); % calculate effective angle of attack at each panel

            % use linear interpolation for C_l and C_d at the effective AoA
            Cl_i = zeros(1, N); % initialise sectional lift coefficient array
            Cd_i = zeros(1, N); % initialise sectional drag coefficient array
            for k = 1:N
                index = find(AoA_data <= AoA_effective(k), 1, 'last'); % find index of last data point below AoA_effective
                % set the index to valid data range
                if isempty(index) || index < 1
                    index = 1; % set the index leftbound 
                elseif index >= length(AoA_data)
                    index = length(AoA_data) - 1; % set the index rightbound
                end
                lin_f = (AoA_effective(k) - AoA_data(index)) / (AoA_data(index+1) - AoA_data(index)); % calculate the linear interpolation factor
                Cl_i(k) = Cl_data(index) + lin_f * (Cl_data(index+1) - Cl_data(index)); % interpolated sectional lift coefficient
                Cd_i(k) = Cd_data(index) + lin_f * (Cd_data(index+1) - Cd_data(index)); % interpolated sectional drag coefficient
            end

            G_new = 0.5 * v_free * c_i .* Cl_i; % calculate new circulation from the Kutta-Joukowski theorem
            G_new = G_i + D * (G_new - G_i); % apply the damping factor

            % check for convergence
            if max(abs(G_new - G_i)) < tol
                converged = true; % change the flag
                G_i = G_new; % set the converged solution
                break % break from the loop
            end

            G_i = G_new; % update circulation for next iteration

        end
        
        % set up a warning call if solution does not converge
        if ~converged
            warning('Solution did not converge for AoA = %.1f deg, N = %d', rad2deg(AoA_array(j)), N)
        end

        l_i = rho * v_free .* G_i; % calculate the lift force per unit span at each panel 
        AoA_induced = AoA_effective - AoA_array(j); % calculate the induced angle at each panel

        Cl_integrand = c_i .* Cl_i .* cos(AoA_induced); % C_L integrand at each panel
        Cd_integrand = c_i .* (Cd_i - Cl_i .* sin(AoA_induced)); % C_D integrand at each panel

        % Simpson's rule integration across all panels
        CL_sum = 0; % initialise C_L integration sum
        CD_sum = 0; % initialise C_D integration sum
        for k = 1:N
            if k == 1 % left tip pannel
                l_left  = 0; % zero lift at left wingtip
                l_mid   = Cl_integrand(k); % value at panel centre
                l_right = (Cl_integrand(k) + Cl_integrand(k+1)) / 2; % average at right boundary
                d_left  = 0; % zero drag integrand at left wingtip
                d_mid   = Cd_integrand(k); % value at panel centre
                d_right = (Cd_integrand(k) + Cd_integrand(k+1)) / 2; % average at right boundary
            elseif k == N % right tip panel
                l_left  = (Cl_integrand(k-1) + Cl_integrand(k)) / 2; % average at left boundary
                l_mid   = Cl_integrand(k); % value at panel centre
                l_right = 0; % zero lift at right wingtip
                d_left  = (Cd_integrand(k-1) + Cd_integrand(k)) / 2; % average at left boundary
                d_mid   = Cd_integrand(k); % value at panel centre
                d_right = 0; % zero drag integrand at right wingtip
            else % panels in between tips
                l_left  = (Cl_integrand(k-1) + Cl_integrand(k)) / 2; % average at left boundary
                l_mid   = Cl_integrand(k); % value at panel centre
                l_right = (Cl_integrand(k) + Cl_integrand(k+1)) / 2; % average at right boundary
                d_left  = (Cd_integrand(k-1) + Cd_integrand(k)) / 2; % average at left boundary
                d_mid   = Cd_integrand(k); % value at panel centre
                d_right = (Cd_integrand(k) + Cd_integrand(k+1)) / 2; % average at right boundary
            end
            CL_sum = CL_sum + (delta_y/6) * (l_left + 4*l_mid + l_right); % apply Simpson's rule for C_L
            CD_sum = CD_sum + (delta_y/6) * (d_left + 4*d_mid + d_right); % apply Simpson's rule for C_D
        end

        CL = CL_sum / S; % normalise C_L by wing area
        CD = CD_sum / S; % normalise C_D by wing area

        % store results for plotting
        CL_results(i, j) = CL;
        CD_results(i, j) = CD;
        
        fprintf('AoA = %.1f deg, N = %d: CL = %.4f, CD = %.4f\n', rad2deg(AoA_array(j)), N, CL, CD) % print results to command window

        fprintf(results_file, '  %-10.1f %-10d %-10.4f %-10.4f\n', rad2deg(AoA_array(j)), N, CL, CD); % append a row of results to the text file

        % plot lift distribution for modes 1 and 3
        if mode == 1 || mode == 3
            if mode == 3 && i == 1
                lift_distribution = figure; % create a figure for mode 3
                hold on; % hold figure since it will contain curves for multiple N values
            elseif mode == 1
                lift_distribution = figure; % create a figure for mode 1
            end

            plot(y_i, l_i,'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', N)); % plot lift distribution
            xlabel('Position along the span / m', 'FontSize', 12); 
            ylabel('Lift per unit span / [N/m]', 'FontSize', 12); 
            title(sprintf('Lift Distribution, AoA = %.1f deg', rad2deg(AoA_array(j))), 'FontSize', 12); 
            legend('FontSize', 12); 
            grid on; 
            xline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off'); % add a vertical axis
            yline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off'); % add a horizontal axis
            if mode == 1
                saveas(lift_distribution, sprintf('lift_distribution_AoA=%.1f_N%d.png', rad2deg(AoA_array(j)), N)); % save the plot for mode 1
            end
        end

    end
end

fclose(results_file); % close the results text file using fclose

if mode == 3
    saveas(lift_distribution, sprintf('lift_distribution_AoA=%.1f_allN.png', rad2deg(AoA_array(1)))); % save the plot for mode 3 including all lines
end

% C_L and C_D vs AoA plots for mode 2
if mode == 2
    AoA_deg = rad2deg(AoA_array); % convert AoA array back to degrees for plotting

    CL_plot = figure; % create a figure for the C_L plot
    plot(AoA_deg, CL_results(1,:), 'LineWidth', 1.5); % plot C_L vs AoA
    xlabel('Angle of Attack / degrees', 'FontSize', 12); 
    ylabel('C_L', 'FontSize', 12); 
    title(sprintf('3D Wing Lift Coefficient, N = %d', N), 'FontSize', 12); 
    grid on;
    xline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
    yline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    saveas(CL_plot, sprintf('CL_vs_AoA_N=%d.png', N)); % save the C_L figure for the given number of panels

    CD_plot = figure; % create a figure for C_D plot
    plot(AoA_deg, CD_results(1,:), 'LineWidth', 1.5); % plot C_D vs AoA
    xlabel('Angle of Attack / degrees', 'FontSize', 12);
    ylabel('C_D', 'FontSize', 12);
    title(sprintf('3D Wing Drag Coefficient, N = %d', N), 'FontSize', 12);
    grid on;
    xline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
    yline(0, 'y', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
    saveas(CD_plot, sprintf('CD_vs_AoA_N=%d.png', N)); % save the C_D figure for the given number of panels
end