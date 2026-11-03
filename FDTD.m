fprintf('=== Convergence Test ===\n');
omega_values_q1 = [25, 50, 100];
N_values_q1 = [40, 80, 160, 320, 640, 1280];
T_q1 = 0.8;
a = 0.25; % Δt = a * Δx

errors_q1 = zeros(length(omega_values_q1), length(N_values_q1));
dx_values_q1 = 1./N_values_q1;


for omega_idx = 1:length(omega_values_q1)
    omega = omega_values_q1(omega_idx);
    
    fprintf('\n=== for ω = %d ===\n', omega);
    
    %% Loop through different grid sizes
    for N_idx = 1:length(N_values_q1)
        N = N_values_q1(N_idx);
        dx = 1/N;
        dy = dx;
        dt = a * dx;
        
       
        M = ceil(T_q1/dt);
        dt = T_q1/M; 
        
        
        % Electric field at integer points
        x_E = 0:dx:1;
        y_E = 0:dy:1;
        
        % Magnetic fields at half-points
        x_Hx = x_E; 
        y_Hx = dy/2:dy:1-dy/2; 
        
        x_Hy = dx/2:dx:1-dx/2; 
        y_Hy = y_E;  
        
        % Initialize fields
        Ez = zeros(length(x_E), length(y_E));      
        Ez_new = zeros(length(x_E), length(y_E));  
        Hx = zeros(length(x_Hx), length(y_Hx));    
        Hy = zeros(length(x_Hy), length(y_Hy));    
        
        % Create mesh for source term
        [X_E, Y_E] = meshgrid(x_E, y_E);
        
        % Initial condition for Ez at t=0
        Ez = 16 * X_E'.^2 .* (X_E'-1).^2 .* Y_E'.^2 .* (Y_E'-1).^2;
        
        % Source term Jz (spatial part)
        Jz = 16*omega * X_E'.^2 .* (X_E'-1).^2 .* Y_E'.^2 .* (Y_E'-1).^2 + ...
             32/omega * ((6*X_E'.^2 - 6*X_E' + 1) .* Y_E'.^2 .* (Y_E'-1).^2 + ...
                        (6*Y_E'.^2 - 6*Y_E' + 1) .* X_E'.^2 .* (X_E'-1).^2);
        
        % Time stepping
        for n = 1:M-1
            t_half = (n + 0.5) * dt;  % Time at n+1/2
            
            % Update Hx using central differences
            for i = 1:length(x_Hx)
                for j = 1:length(y_Hx)
                  
                    
                    y_idx_j = j;  % Corresponds to Ez at y = (j-1)*dy
                    y_idx_jp1 = j+1;  % Corresponds to Ez at y = j*dy
                    
                    if y_idx_jp1 <= length(y_E)
                        Ez_bottom = Ez(i, y_idx_jp1);
                    else
                        Ez_bottom = 0;  % Isolated boundary condition
                    end
                    
                    Ez_top = Ez(i, y_idx_j);
                    
                    % Central difference for ∂Ez/∂y
                    dEz_dy = (Ez_bottom - Ez_top) / dy;
                    
                    % Update Hx
                    Hx(i, j) = Hx(i, j) - dt * dEz_dy;
                end
            end
            
            % Update Hy  using central differences
            for i = 1:length(x_Hy)
                for j = 1:length(y_Hy)
                   
                    
                    x_idx_i = i;  % Corresponds to Ez at x = (i-1)*dx
                    x_idx_ip1 = i+1;  % Corresponds to Ez at x = i*dx
                    
                    if x_idx_ip1 <= length(x_E)
                        Ez_right = Ez(x_idx_ip1, j);
                    else
                        Ez_right = 0;  % Beyond boundary
                    end
                    
                    Ez_left = Ez(x_idx_i, j);
                    
                    % Central difference for ∂Ez/∂x
                    dEz_dx = (Ez_right - Ez_left) / dx;
                    
                    % Update Hy
                    Hy(i, j) = Hy(i, j) + dt * dEz_dx;
                end
            end
            
            % Update Ez 
            for i = 2:length(x_E)-1  % Interior points only
                for j = 2:length(y_E)-1
       
                    if j-1 >= 1
                        Hx_top = Hx(i, j-1);
                    else
                        Hx_top = 0;  % Boundary
                    end
                    
                    if j <= length(y_Hx)
                        Hx_bottom = Hx(i, j);
                    else
                        Hx_bottom = 0;  % Boundary
                    end
                    
                    dHx_dy = (Hx_bottom - Hx_top) / dy;
                    
                 
                    if i-1 >= 1
                        Hy_left = Hy(i-1, j);
                    else
                        Hy_left = 0;  % Boundary
                    end
                    
                    if i <= length(x_Hy)
                        Hy_right = Hy(i, j);
                    else
                        Hy_right = 0;  % Boundary
                    end
                    
                    dHy_dx = (Hy_right - Hy_left) / dx;
                    
                    % Update Ez
                    curl_H = dHy_dx - dHx_dy;
                    Ez_new(i, j) = Ez(i, j) + dt * (curl_H - sin(omega * t_half) * Jz(i, j));
                end
            end
            
       
            Ez_new(1, :) = 0;     
            Ez_new(end, :) = 0;    
            Ez_new(:, 1) = 0;     
            Ez_new(:, end) = 0;   
            
            % Prepare for next time step
            Ez = Ez_new;
        end
        
        % Compute exact solution
        Ez_exact = 16 * X_E'.^2 .* (X_E'-1).^2 .* Y_E'.^2 .* (Y_E'-1).^2 * cos(omega * T_q1);
        
        % Compute maximum error
        Ez_interior = Ez(2:end-1, 2:end-1);
        Ez_exact_interior = Ez_exact(2:end-1, 2:end-1);
        max_error = max(abs(Ez_interior(:) - Ez_exact_interior(:)));
        errors_q1(omega_idx, N_idx) = max_error;
        
        fprintf('N = %4d, dx = %.6f, dt = %.6f, Error = %.6e\n', ...
                N, dx, dt, max_error);
    end
end


figure('Position', [100, 100, 1200, 400], 'Name', 'Convergence Test');

markers = {'o-', 's-', '^-'};
colors = {'b', 'r', 'g'};

for omega_idx = 1:length(omega_values_q1)
    subplot(1, 3, omega_idx);
    loglog(N_values_q1, errors_q1(omega_idx, :), markers{omega_idx}, ...
           'Color', colors{omega_idx}, 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    
    ref_order1 = errors_q1(omega_idx, 1) * (N_values_q1(1)./N_values_q1);
    ref_order2 = errors_q1(omega_idx, 1) * (N_values_q1(1)./N_values_q1).^2;
    
    loglog(N_values_q1, ref_order1, 'k--', 'LineWidth', 1.5);
    loglog(N_values_q1, ref_order2, 'k:', 'LineWidth', 1.5);
    
    xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Maximum Error', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$\\omega = %d$', omega_values_q1(omega_idx)), ...
          'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    
    
    rates = log(errors_q1(omega_idx, 1:end-1)./errors_q1(omega_idx, 2:end)) / log(2);
    fprintf('\nω = %d convergence rates:\n', omega_values_q1(omega_idx));
    for i = 1:length(rates)
        fprintf('  N=%d to N=%d: order = %.3f\n', ...
                N_values_q1(i), N_values_q1(i+1), rates(i));
    end
    
    legend('Yee scheme', '1st order', '2nd order', ...
           'Location', 'southwest', 'Interpreter', 'latex');
end

sgtitle('Convergence of Yee Scheme for 2D TM Model', ...
        'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');


fprintf('\n\n=== Wave Propagation with Gaussian Source ===\n');


omega_values_q2 = [10, 50, 100];
N_values_q2 = [40, 80, 160];
T_q2 = 10; 
a_q2 = 0.25; 
plot_times = [2.5, 5, 7.5, 10];

x_min = -1; x_max = 1;
y_min = -1; y_max = 1;


for omega_idx = 1:length(omega_values_q2)
    omega = omega_values_q2(omega_idx);
    
    fprintf('\n=== for ω = %d ===\n', omega);
    
   
    figure('Position', [50, 50, 1400, 900], 'Name', sprintf('ω = %d', omega));
    
   
    sigma = max(36, omega^2);
    
    for N_idx = 1:length(N_values_q2)
        N = N_values_q2(N_idx);
        
        % Spatial step
        dx = (x_max - x_min)/N;
        dy = dx;
        dt = a_q2 * dx;
        
        % Number of time steps
        M = ceil(T_q2/dt);
        dt = T_q2/M;
        
        % Grid points
        x = x_min:dx:x_max;
        y = y_min:dy:y_max;
        [X, Y] = meshgrid(x, y);
        
        % Initialize fields
        Ez = zeros(length(x), length(y));
        Ez_new = zeros(length(x), length(y));
        Hx = zeros(length(x), length(y)-1);  % Staggered in y
        Hy = zeros(length(x)-1, length(y));  % Staggered in x
        
        
        Jz_spatial = -omega * exp(-sigma * ((X' - 0.01).^2 + (Y' - 0.015).^2));
        
        % Store snapshots
        snapshots = cell(length(plot_times), 1);
        time_steps_to_save = round(plot_times / dt);
        
        % Time stepping
        for n = 0:M-1
            t = n * dt;
            t_half = (n + 0.5) * dt;
            
          
            [is_save, idx] = ismember(n, time_steps_to_save);
            if is_save
                snapshots{idx} = Ez;
            end
            
            % Update Hx
            for i = 1:length(x)
                for j = 1:length(y)-1
                    % Central difference for ∂Ez/∂y
                    dEz_dy = (Ez(i, j+1) - Ez(i, j)) / dy;
                    Hx(i, j) = Hx(i, j) - dt * dEz_dy;
                end
            end
            
            % Update Hy
            for i = 1:length(x)-1
                for j = 1:length(y)
                    % Central difference for ∂Ez/∂x
                    dEz_dx = (Ez(i+1, j) - Ez(i, j)) / dx;
                    Hy(i, j) = Hy(i, j) + dt * dEz_dx;
                end
            end
            
            % Update Ez
            for i = 2:length(x)-1
                for j = 2:length(y)-1
                    % Calculate curl of H = ∂Hy/∂x - ∂Hx/∂y
                    
                    % ∂Hy/∂x: (Hy(i,j) - Hy(i-1,j))/dx
                    if i-1 >= 1 && i <= size(Hy,1)
                        dHy_dx = (Hy(i, j) - Hy(i-1, j)) / dx;
                    else
                        dHy_dx = 0;
                    end
                    
                    % ∂Hx/∂y: (Hx(i,j) - Hx(i,j-1))/dy
                    if j-1 >= 1 && j <= size(Hx,2)
                        dHx_dy = (Hx(i, j) - Hx(i, j-1)) / dy;
                    else
                        dHx_dy = 0;
                    end
                    
                    curl_H = dHy_dx - dHx_dy;
                    Ez_new(i, j) = Ez(i, j) + dt * (curl_H - sin(omega * t_half) * Jz_spatial(i, j));
                end
            end
            
            
            Ez_new(1, :) = 0;      
            Ez_new(end, :) = 0;    
            Ez_new(:, 1) = 0;      
            Ez_new(:, end) = 0;    
            Ez = Ez_new;
        end
        
      
        if isempty(snapshots{end})
            snapshots{end} = Ez;
        end
        
       
        for t_idx = 1:length(plot_times)
            subplot(length(N_values_q2), length(plot_times), (N_idx-1)*length(plot_times) + t_idx);
            
            if ~isempty(snapshots{t_idx})
                imagesc(x, y, snapshots{t_idx}');
                axis xy; axis equal; axis tight;
                colorbar;
                colormap(jet);
                
                xlabel('x');
                ylabel('y');
                title(sprintf('N=%d, t=%.1f', N, plot_times(t_idx)));
                clim([-1, 1]);  
            end
        end
    end
    
    sgtitle(sprintf('E_z Field Evolution for ω = %d', omega), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    fprintf('  Completed ω = %d\n', omega);
end