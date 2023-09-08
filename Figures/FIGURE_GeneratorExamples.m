function FIGURE_GeneratorExamples
% This function generates the figure showing different patterns that the
% generator is capable of producing

% Define the dimensions of the images to plot
Nx = 500;
Ny = 500;

% Define the parameter values used for the three different rows of the
% figure
fibre_sep = [0.6, 0.25, 1.5];
feature_size = [0.08, 0.8, 2];
roughness = [0.1, 0.4, 0.9];
patch_size = [2, 8, 3];
anisotropy = [1, 5, 1];
direction_vec = [45, 15, 75];

% Specify the two different combinations of fibreness and patchiness that
% define the patterns in the two rightmost columns
patchiness = [0.05, 0.5];
fibreness = [0.4, 0.05];

% Define the plotting sizes (all as fractions of one)
margin = 0.025;
xgap = -0.07;
xsep = 0.025;
ygap = 0.02;
leftTextSpace = 0.1;
titleSpace = 0.1;

% Define density
density = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch up a level to the main folder to run code
cd ..

% Load in seed data
load('fibro_seedinfo.mat','permute_tables','offset_tables');
% Just use the first few seeds
Ps = {permute_tables{1}, permute_tables{2}, permute_tables{3}};
offsets = {offset_tables{1}, offset_tables{2}, offset_tables{3}};

% Define the mesh
mesh = buildMesh(Nx,Ny,1/136);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the number of rows and columns from the parameter vectors
n_rows = length(fibre_sep);
n_cols = length(patchiness) + 3;   % +3 for visualisation of the 3 fields

% Calculate the derived plotting sizes
dx = (1 - 2*margin - (n_cols-1)*xgap - leftTextSpace - xsep) / n_cols;
dy = (1 - 2*margin - (n_rows-1)*ygap - titleSpace) / n_rows;

% Define final plot titles
field_texts = {'Base Noise, $\mathcal{N}_b$', 'Density Variation, $\mathcal{N}_d$', 'Fibre Selection, $\mathcal{F}$'};
pattern_titles = {'$$f = 0.4,\,\,\,\,\,d = 0.05$$', '$$f = 0.05,\,\,\,\,\,d = 0.5$$'};

% Convert the direction vector to radians
direction = direction_vec * pi / 180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a 'fibrosis' colormap for the resulting patterns
fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]];

% Load in the 'viridis' colormap for the noisefields
load('extra_colormaps.mat','viridis');

% Initialise a new figure, fullscreen
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Create a series of axes objects that will be used for plotting
counter = 0;
for k = 1:n_rows           % rows
    for j = 1:n_cols       % columns

        % Increment axes counter
        counter = counter + 1;
        
        % Create axes object
        xpos = margin + leftTextSpace + (j-1)*(xgap + dx) + (j > 3)*xsep;
        ypos = margin + (n_rows - k)*(ygap + dy);
        ax{k,j} = axes('Position', [xpos, ypos, dx, dy]);
        
    end
end

% Now actually do the plotting
for k = 1:n_rows    
    
    % Create the noise fields, and the patterns for this set of parameters
    for m = 1:length(fibreness)
        params = [fibreness(m), fibre_sep(k), patchiness(m), feature_size(k), roughness(k), patch_size(k), anisotropy(k), direction(k)];
        [patterns{m}, O_b, O_d, F] = createFibroPattern(mesh, density, params, Ps{k}, offsets{k});
    end
    
    for j = 1:n_cols
        
        if j <= 3    % Plotting of noise fields
            switch j
                case 1    % Plot base fibrotic field
                    imagesc(ax{k,j}, O_b); axis(ax{k,j}, 'equal', 'off');
                    colormap(ax{k,j}, viridis);
                case 2    % Plot fibre field
                    imagesc(ax{k,j}, O_d); axis(ax{k,j}, 'equal', 'off');
                    colormap(ax{k,j}, viridis);
                case 3    % Plot density variation field
                    imagesc(ax{k,j}, F); axis(ax{k,j}, 'equal', 'off');
                    colormap(ax{k,j}, viridis);
            end
            
            if k == 1   % Titles on top row
                text(ax{k,j}, Nx/2,-Ny/10, field_texts{j}, 'FontSize', 20, 'HorizontalAlignment', 'Center', 'Interpreter', 'latex');
            end
            
        else         % Plotting of final patterns
            
            imagesc(ax{k,j}, patterns{j-3}); axis(ax{k,j}, 'equal', 'off');
            colormap(ax{k,j},fibroclr);
            
            if k == 1   % Titles on top row
                text(ax{k,j}, Nx/2,-3*Ny/40, pattern_titles{j-3}, 'FontSize', 20, 'HorizontalAlignment', 'Center', 'Interpreter', 'latex');
            end
            
        end 
    end
    
    % Specify the variable values on the left of the figure
    text(ax{k,1}, -Nx/30, Ny/2, ['$\begin{array}{lll}l_b = ', num2str(feature_size(k)),' & & r = ', num2str(anisotropy(k)), '\\ \gamma = ', num2str(roughness(k)), '& & l_d = ',num2str(patch_size(k)),'\\ L = ',num2str(fibre_sep(k)), '& & \phi = ', num2str(direction_vec(k)),'^{\circ} \end{array}$'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'Fontsize', 20);
    
    
end

% Manually place mega-titles
text(ax{1,2}, Nx/2, -3*Ny/10, 'Component Noise Fields', 'FontSize', 36, 'HorizontalAlignment', 'Center');
text(ax{1,4}, Nx * 1.05, -3*Ny/10, 'Resulting Patterns', 'FontSize', 36, 'HorizontalAlignment', 'Center');
           
% Return to figures folder
cd Figures