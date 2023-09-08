function FIGURE_VisualiseEllipseMetrics
% This function generates the figure showing the four histological sections
% that demonstrate the different classifications of microfibrosis

% Define the plotting sizes (all as fractions of one)
margin = 0.075;
xgap = 0.02;
ygap = 0.04;
titleSpace = 0.05;
leftWidth = 0;

x_size = 0.6;
y_size = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch up a level to the main folder to run code
cd ..

% Load the representative histological sections
load('histo_patterns.mat','patterns');

% Read out the dimensions of these patterns
[Ny, Nx] = size(patterns{1});

% Calculate the metrics associated with each pattern
for k = 1:length(patterns)
    metrics{k} = ellipseMetrics(patterns{k}, [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of rows and columns from the parameter vectors
n_rows = 2;                  % Two rows - the patterns, and the ellipses
n_cols = length(patterns);   % One column per pattern

% Calculate the derived plotting sizes
dx = (1 - 2*margin - (n_cols-1)*xgap - leftWidth) / n_cols;
dy = (1 - 2*margin - titleSpace);

% Define final plot titles
fibro_classes = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise a new figure
figure('units', 'normalized', 'OuterPosition', [0.5-x_size/2 0.5-y_size/2 0.5+x_size/2 0.5+y_size/2]);

% Create a series of axes objects that will be used for plotting
counter = 0;
for j = 1:n_cols       
    
    % Increment axes counter
    counter = counter + 1;
    
    % Create axes object
    xpos = margin + leftWidth + (j-1)*(xgap + dx);
    ypos = 1 - margin - dy;
    ax{j} = axes('Position', [xpos, ypos, dx, dy]);
    
end

% Now actually do the plotting - loop over columns
for j = 1:n_cols
        
    % Plot the ellipses corresponding to each histological pattern
    fig = gcf;
    fig.CurrentAxes = ax{j};
    plotEllipses(metrics{j});
    
end

% Re-set the axes for the ellipses to be consistent across all four (to
% allow for a consistent visual comparison)

% Find the maximum ellipse size - loop over cell elements and combine all
% into one big long vector to check for maximum
check_vec = [];
for j = 1:n_cols
    check_vec = [check_vec; metrics{j}(2:3:end)'];       % 2nd metric for each ellipse is the major axis
end
max_size = max(check_vec);

for j = 1:n_cols
    
    % Set axis to the maximum size, with a little buffer
    axis(ax{j}, 0.55*[-max_size max_size -max_size max_size]);
    axis(ax{j}, 'square');
    set(ax{j},'XTick',[],'YTick',[],'LineWidth',2)
    box(ax{j},'on');
    
    % Also add titles here, because axis dimensions are now known
    text(ax{j}, 0, max_size * 0.65, fibro_classes{j}, 'HorizontalAlignment', 'Center', 'FontSize', 24);
    
end

% Return to figures folder
cd Figures