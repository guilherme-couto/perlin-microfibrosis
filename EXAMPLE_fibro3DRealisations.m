function mesh3D = EXAMPLE_fibro3DRealisations(seed_no)
% This function creates an example pattern of the four types of fibrosis,
% using the single parameter sets found by calibrating the generator to the
% two-dimensional histological images.
%
% Usage:
%     field = EXAMPLE_fibro3DRealisations(seed_no)
%
% INPUTS:
% 
% seed_no: Seed for randomness used in generating a specific noisefield
%          realisation (7 was used in all cases in the paper involving 
%          selecting a single seed)
%
% Additional properties of the generation can be edited at the top of the
% .m file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define dimensions (millimetres)
Lx = 1.5;
Ly = 1.5;
Lz = 0.4;                     % third, "transmural" direction

% Define the side length of the 'voxels' (millimetres)
dx = 1/136;       

% Specify the modifier to anisotropy in the third dimension
K = 1;                         % no modification

% Specify the density correction function representing shift to 3D
% (see paper supplement for justification)
density_fun = @(rho) rho * 2/3;

% Define the pattern names in order (as they appear in stored generator
% calibrations file)
pattern_names = {'interstitial', 'compact', 'diffuse', 'patchy'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% MESH SETUP

% Determine the number of points in each direction
Nx = floor( Lx / dx ) + 1;
Ny = floor( Ly / dx ) + 1;
Nz = floor( Lz / dx ) + 1;

% Create a three-dimensional mesh object
[X,Y,Z] = meshgrid(0:dx:Lx, 0:dx:Ly, 0:dx:Lz);
points = [X(:) - mean(X(:)), Y(:) - mean(Y(:)), Z(:) - mean(Z(:))]';
mesh3D.N_points = [Nx, Ny, Nz];
mesh3D.points = points';


%%% GENERATOR SETUP

% Load in the permutation tables
load('fibro_seedinfo3d.mat','permute_tables','offset_tables');

% Load in the generator calibrations learned from 2D histology
load('param_modes.mat','param_modes','densities');

% Convert parameters as needed in moving from 2D to 3D
feature_size3 = param_modes(:,4) .* param_modes(:,7).^(-1/6) * K^(-1/3);
anisotropy_y = param_modes(:,7);
anisotropy_z = param_modes(:,7) * K;

% Set all angles to zero for simplicity
phi = zeros(4,1);
theta = zeros(4,1);

% Convert from the 2D learned parameters to a 3D equivalent
params3D = [param_modes(:,1), param_modes(:,2), param_modes(:,2), param_modes(:,3), feature_size3, param_modes(:,5), param_modes(:,6), anisotropy_y, anisotropy_z, phi, theta];

% Grab out information corresponding to input seed - two different sets of
% offset data are needed for three dimensions
Ps = permute_tables{seed_no};
offsets = offset_tables{seed_no};
offsets2 = offset_tables{seed_no+1};


%%% PATTERN CREATION AND SAVING

% Loop over the different types of fibrosis
for k = 1:size(param_modes,1)
       
    % Apply the density correction for three dimensions to the density
    % learned for this pattern in 2D
    density = density_fun( densities(k) );
    % Use the 3D version of the generator to get a binary pattern
    presence = createFibroPattern3D(mesh3D, density, params3D(k,:), Ps, offsets, offsets2);
    % Convert this file into a VTK using pattern name as filename
    writeVTKDataNew(mesh3D, presence, pattern_names{k});
    
end
