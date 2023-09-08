function FIGURE_PatternMatchParams
% This function plots the figure that shows the parameter values selected
% for the four histological patterns

% Specify the populations to load
particle_filenames = {'interstitial2000_full', 'compact2000_full_nofibres', 'diffuse2000_full_nofibres', 'patchy2000_full'};
pattern_names = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};

% Variable ordering (for plotting)
variable_ordering = [4, 5, 7, 8, 3, 6, 1, 2];
variable_names = {'Fibreness', 'Fibre Separation', 'Density Var.', 'Feature Size', 'Roughness', 'Density Var. Scale', 'log Anisotropy', 'Direction'};

% Variables present
present_vars = { [ 1, 2, 3, 4, 5, 6, 7, 8]; [3, 4, 5, 6, 7, 8]; [3, 4, 5, 6, 7, 8]; [1, 2, 3, 4, 5, 6, 7, 8]};

% Variable ranges
theta_mins = [ 0, 0.3, 0, 0.01, 0, 1, 1, -pi/2 ];
theta_maxs = [ 0.4, 2, 0.5, 2, 0.99, 8, 50, pi/2 ];
scale_param = logical([ 0, 0, 0, 0, 0, 0, 1, 0]);

% Define axis positions, etc.
x_margin = 0.025;
y_margin = 0.1;
xgap = 0.025;
ygap = 0.05;
leftspace = 0.075;
topspace = -0.075;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First change folder back to the base folder
cd ..

% Read out number of rows of results to plot
n_patterns = length(particle_filenames);

% Load all particle sets into data initially
for j = 1:n_patterns
    load(['Results\',particle_filenames{j},'.mat'],'particles');
    particle_sets{j} = particles;
end
n_vars = length(theta_mins);

% Transform the scale parameters' ranges by logarithm
theta_mins(scale_param) = log(theta_mins(scale_param));
theta_maxs(scale_param) = log(theta_maxs(scale_param));

% Create the linspaces used for plotting the density estimates
n_points = 500;
theta_vals = zeros(n_vars, n_points);
for i = 1:n_vars
   theta_vals(i,:) = linspace(theta_mins(i), theta_maxs(i), n_points);
end

% Convert theta range for angles to degrees
theta_mins(8) = theta_mins(8) * 180/pi;
theta_maxs(8) = theta_maxs(8) * 180/pi;
theta_vals(8,:) = theta_vals(8,:) * 180/pi;

% Calculate the axis sizes
dx = (1 - 2*x_margin - leftspace - (n_vars-1) * xgap) / n_vars;
dy = (1 - 2*y_margin - topspace - (n_patterns-1) * ygap) / n_patterns;

% Initialise figure
figure('Units','normalized','OuterPosition',[0 0 1 1]);

% Set up the axes for each individual plot
for j = 1:n_patterns
    for i = 1:length(present_vars{j})
        ax{j,i} = axes('Position',[x_margin + leftspace + (i-1) * (xgap+dx), 1 - y_margin - topspace - (j) * (ygap+dy), dx, dy]);
    end
end

% Loop over each variable and plot it, but limited only to the variables
% present for this particle set
for j = 1:n_patterns
    for i = 1:n_vars
        
        % Read out the variable to plot in this location (from ordering)
        plot_var = variable_ordering(i);
        
        % Check if this variable is present for this dataset
        loc = find(present_vars{j} == plot_var);
        if loc
            
            % If this is direction, convert to degrees
            if plot_var == 8
                particle_sets{j}.thetas(:,loc) = particle_sets{j}.thetas(:,loc) * 180/pi;
            end
            
            % Create the kernel-estimated density
            [p, theta] = ksdensity(particle_sets{j}.thetas(:,loc), theta_vals(plot_var,:));
            
            % Plot this on the required axes
            plot( ax{j,i}, theta, p, 'LineWidth', 2)
            
            % Set axis limits
            axis( ax{j,i}, [theta_mins(plot_var) theta_maxs(plot_var) 0 1.1*max(p)]);
            
            % Increase fontsize of ticks
            set(ax{j,i}, 'FontSize', 14);
            
        end
        
        % If this is the first column, add text on the left (pattern name)
        if i == 1
           text( ax{j,i}, theta_mins(plot_var) - 0.2 * (theta_maxs(plot_var) - theta_mins(plot_var)), 0.55 * max(p), pattern_names{j}, 'HorizontalAlignment', 'Right', 'FontSize', 20);
        end
        
        % If this is the first row, add text above (variable name)
        if j == 1
            text( ax{j,i}, theta_mins(plot_var) + 0.5 * (theta_maxs(plot_var) - theta_mins(plot_var)), 1.3*max(p), variable_names{plot_var}, 'HorizontalAlignment', 'Center', 'FontSize', 20);
        end
        
    end
    
end

% Return to figures folder
cd Figures

end

