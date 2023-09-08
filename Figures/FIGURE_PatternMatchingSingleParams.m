function FIGURE_PatternMatchingSingleParams
% This function plots some patterns generated using the proposed generator
% tunings in the paper

% Define the plotting sizes (all as fractions of one)
margin = 0.025;
xgap = -0.25;
xsep = 0.03;
ygap = 0.02;
leftTextSpace = 0.025;
titleSpace = 0.05;

% Specify the number of generated patterns to show
n_show = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch up a level to the main folder to run code
cd ..

% Define the fibrosis classification names
fibro_classes = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};

% Calculate the number of rows and columns from the parameter vectors
n_rows = length(fibro_classes);
n_cols = n_show + 1;                  % +1 for the target pattern

% Calculate the derived plotting sizes
dx = (1 - 2*margin - (n_cols-1)*xgap - leftTextSpace - xsep) / n_cols;
dy = (1 - 2*margin - (n_rows-1)*ygap - titleSpace) / n_rows;

% Define a 'fibrosis' colormap
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Initialise a new figure, fullscreen
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Create a series of axes objects that will be used for plotting
counter = 0;
for k = 1:n_rows           % rows
    for j = 1:n_cols       % columns

        % Increment axes counter
        counter = counter + 1;
        
        % Create axes object
        xpos = margin + leftTextSpace + (j-1)*(xgap + dx) + (j>1) * xsep;
        ypos = margin + (n_rows - k)*(ygap + dy);
        ax{k,j} = axes('Position', [xpos, ypos, dx, dy]);
        
    end
end

% Load the patterns from the histological images
load('histo_patterns.mat','patterns');

% Read out the dimensions of the patterns
[Ny, Nx] = size(patterns{1});

% Load in the generator tunings
load('param_modes.mat','param_modes');

% Load in the seed information
load('fibro_seedinfo.mat');

% Create the mesh
mesh = buildMesh(250, 400, 1/136);

% Loop over the classes of microfibrosis, plotting the original
% histological image, and then the generated patterns
for k = 1:length(fibro_classes)
   
    % First, plot the histological section
    imagesc(ax{k,1}, patterns{k}); axis(ax{k,1}, 'equal', 'off');
    colormap(ax{k,1}, fibroclr);
    
    % Calculate its density (for generating the following patterns)
    density = sum(sum(patterns{k}))/numel(patterns{k});
    
    for j = 1:n_show
        
        % Generate the pattern for the current class and seed number
        if any(isnan(param_modes(k,:)))
            [presence, ~, ~] = createFibroPatternNoFibres(mesh, density, param_modes(k,~isnan(param_modes(k,:))), permute_tables{j}, offset_tables{j});
        else
            [presence, ~, ~, ~] = createFibroPattern(mesh, density, param_modes(k,:), permute_tables{j}, offset_tables{j});
        end
        
        % Plot this pattern
        imagesc(ax{k,j+1}, presence); axis(ax{k,j+1}, 'equal', 'off');
        colormap(ax{k,j+1}, fibroclr);
        
    end
        
end

% Add column titles
text(ax{1,1}, Nx/2, -Ny/6, 'Histological', 'HorizontalAlignment', 'Center', 'FontSize', 24);
text(ax{1,1+round(n_show/2)}, Nx*1.2, -Ny/6, 'Generated', 'HorizontalAlignment', 'Center', 'FontSize', 24);

% Add row titles
for k = 1:length(fibro_classes)
    text(ax{k,1}, -Nx/5, Ny/2, fibro_classes{k}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle', 'FontSize', 24);
end
           
% Return to figures folder
cd Figures