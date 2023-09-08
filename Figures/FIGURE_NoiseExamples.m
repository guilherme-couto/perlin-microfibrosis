function FIGURE_NoiseExamples
% This function visualises different noisefields, to use as examples in the
% paper


% Define the plotting sizes (all as fractions of one)
margin = 0.01;
xgap = -0.25;
rowhead = 0.08;

% Specify the number of rows and columns in which to plot the examples
n_rows = 2;
n_cols = 3;

% Define the example names
names = {'White Noise','Single Octave Perlin Noise','Multiple Octave Perlin Noise','Anisotropic Octave Noise','Additive Combination','Gaussian Random Field'};
identifiers = {'i', 'ii', 'iii', 'iv', 'v', 'vi'};

% Specify the size of the example patterns to generate
Nx = 325;
Ny = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch up a level to the main folder to run code
cd ..

% Create a mesh of points with spacing one (units are just pixels)
mesh = buildMesh(Nx,Ny,1);

% Load up the seed information
load('fibro_seedinfo.mat', 'permute_tables', 'offset_tables');

% Calculate the derived plotting sizes:
dx = (1 - 2*margin - (n_cols - 1) * xgap ) / n_cols;
dy = (1 - 2*margin - n_rows * rowhead ) / n_rows;

% Initialise the figure and axes for plotting
figure('units', 'normalized', 'OuterPosition', [0 0 1 1]);
counter = 0;
for j = 1:n_rows
    for i = 1:n_cols
        
        % Increment axis counter
        counter = counter + 1;
        
        % Specify the axis positions
        xpos = margin + (i-1) * (dx + xgap);
        ypos = 1 - margin - j * (dy + rowhead);
        
        % Create an axis object here
        ax{counter} = axes('Position', [xpos, ypos, dx, dy]);
        
    end
end

% Now create the actual example patterns. This is not a for loop, because
% each pattern is generated in its own unique way

% EXAMPLE ONE - White Noise
patterns{1} = rand(Ny,Nx);

% EXAMPLE TWO - Perlin Noise
feature_size = 48;
noisefield = Octave2D( mesh.points' / feature_size, 1, 0.5, permute_tables{2}, offset_tables{2} );
patterns{2} = reshape(noisefield, Nx, Ny)';

% EXAMPLE THREE - Octave Noise
feature_size = 48;
noisefield = Octave2D( mesh.points' / feature_size, 4, 0.5, permute_tables{2}, offset_tables{2} );
patterns{3} = reshape(noisefield, Nx, Ny)';

% EXAMPLE FOUR - Anisotropic Octave Noise
feature_size = 48;
theta = pi/3;                % Rotation angle
r = 3.5;                       % Anisotropy ratio
R = [[cos(theta), sin(theta)]; [-sin(theta), cos(theta)]];
A = diag([1/sqrt(r), sqrt(r)]);
eval_points = mesh.points * R' * A / feature_size;
noisefield = Octave2D( eval_points', 4, 0.5, permute_tables{2}, offset_tables{2} );
patterns{4} = reshape(noisefield, Nx, Ny)';

% EXAMPLE FIVE - Density Variation by Additive Combination
feature_size = 48;
feature_size2 = 400;
theta = pi/3;                % Rotation angle
r = 3.5;                       % Anisotropy ratio
R = [[cos(theta), sin(theta)]; [-sin(theta), cos(theta)]];
A = diag([1/sqrt(r), sqrt(r)]);
eval_points = mesh.points * R' * A / feature_size;
noisefield = Octave2D( eval_points', 4, 0.5, permute_tables{2}, offset_tables{2} );
noisefield2 = Octave2D( mesh.points' / feature_size2, 4, 0.5, permute_tables{10}, offset_tables{10} );
patterns{5} = reshape(0.6 * noisefield + 0.4*noisefield2, Nx, Ny)';


% EXAMPLE SIX - Gaussian Random Field

% (Removed for Github version, as this uses code form the cited work and
% this avoids any licensing issues
set(ax{6},'Visible',false);

% Load in the colormap
load('extra_colormaps.mat','viridis');

% Now display all patterns
for k = 1:length(patterns)
    
    % Visualise pattern
    imagesc(ax{k}, flipud(patterns{k}));
    
    % Set axis properties and colormap
    axis(ax{k}, 'off', 'image');
    colormap(ax{k}, viridis);
    
    % Add the name of this pattern above it
    text(ax{k}, Nx/2, -Ny/15, ['(', identifiers{k}, ') ', names{k}], 'FontSize', 24, 'HorizontalAlignment', 'Center');
    
end

% Return to figures folder
cd Figures

end