function FIGURE_VisualiseHistology
% This function generates the figure showing the four histological sections
% that demonstrate the different classifications of microfibrosis

% Define the plotting sizes (all as fractions of one)
margin = 0.075;
xgap = 0.01;
ygap = 0.04;
titleSpace = 0.05;
leftWidth = 0;

x_size = 0.6;
y_size = 0.4;

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
n_cols = length(patterns);   % One column per pattern

% Calculate the derived plotting sizes
dx = (1 - 2*margin - (n_cols-1)*xgap - leftWidth) / n_cols;
dy = (1 - 2*margin - titleSpace);

% Define final plot titles
fibro_classes = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Initialise a new figure, fullscreen
figure('units', 'normalized', 'OuterPosition', [0.5-x_size/2 0.5-y_size/2 0.5+x_size/2 0.5+y_size/2]);

% Create a series of axes objects that will be used for plotting
counter = 0;
for j = 1:n_cols       % columns
    
    % Create axes object
    xpos = margin + leftWidth + (j-1)*(xgap + dx);
    ypos = 1 - margin - dy;
    ax{j} = axes('Position', [xpos, ypos, dx, dy]);
    
end

% Now actually do the plotting - loop over columns
for j = 1:n_cols
    
    % Plot the binary versions of the histological images
    imagesc(ax{j}, patterns{j}); axis(ax{1,j}, 'equal', 'off');
    colormap(ax{j}, fibroclr);
    
    % Add titles to each column
    text(ax{j}, Nx/2, -3*Ny/40, fibro_classes{j}, 'HorizontalAlignment', 'Center', 'FontSize', 24);
        
end

% Return to figures folder
cd Figures