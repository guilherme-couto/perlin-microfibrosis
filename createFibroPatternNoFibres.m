function [presence, O_b, O_d] = createFibroPatternNoFibres(mesh, density, params, Ps, offsets)
% This function takes an input mesh of points and generates a binary
% pattern of presence/absence of fibrotic deposition at those points, in 
% effect forming an image. This version neglects the fibre-selecting field.
% Could also be achieved by using the normal creation function and setting
% fibreness to value zero.
%
% Usage:    presence = CreateFibroPatternNoFibres(mesh, density, params, Ps, offsets)
%
% INPUTS:   mesh:      a mesh object (use createMesh.m to create)
%           density:   the proportion of points to be occupied by fibrosis
%           params:    a set of parameters defining the properties of the noise pattern (see below)
%           Ps:        an m x n matrix, which is m rows of the numbers 0:n-1 arranged in random order (permutation tables for random assignment of vectors in Perlin noise)
%           offsets:   an m x 2 matrix that specifies grid offsets for each octave in octave noise
%
%           * m is the maximum number of octaves that will be requested by
%             this function, and n is the number of Perlin vectors - here
%             256 vectors are used
%           * inputs "Ps" and "offsets" can be generated using function 
%             generateSeeddata.m
%
% OUTPUTS:  presence:   a presence/absence map of fibrosis
%           (O_b):      the Perlin noise field for fibrosis (optional)
%           (O_d):      the Perlin noise field for density variation (optional)
%
% PARAMS:   Parameters are provided as a single vector, specified as:
%           [ patchiness, feature_size, roughness, patch_size, alignment_ratio, direction ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INITIAL SETUP

% Read out paramaters from vector
params = num2cell(params);
[patchiness, feature_size, roughness, patch_size, fibre_alignment, direction] = deal(params{:});

% Create a rotated set of points for the application of anisotropy and
% creation of fibre-aligned pattern. Stored as two rows for ease of matrix
% transforms and input into C++ functions
R_points = [ [ cos(direction) sin(direction) ]; [-sin(direction), cos(direction)] ] * mesh.points';

% Create a new permutation table from the provided by applying it to itself
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k,Ps(k,:)+1);
end


%%% MAIN FIBROSIS PATTERNING FIELD

% Transform points according to input parameters, then call Octave2D
P_f_points = [ R_points(1,:) / sqrt(fibre_alignment); R_points(2,:) * sqrt(fibre_alignment) ];
O_b = Octave2D( P_f_points / feature_size, 4, roughness, Ps, offsets);


%%% LARGE-SCALE DENSITY VARIATION

% Use Octave2D with scaling of point co-ords to attain desired patch_size
O_d = Octave2D(mesh.points' / patch_size, 3, 0.5, Ps2, offsets);


%%% FIBROSIS PATTERN CREATION

% Final noisefield is the weighted combination of the component fields
noise = O_b + patchiness * O_d;

% Convert the final noisefield into a binary pattern
presence = thresholdPattern(noise, density);

% Shift from a list of points back into a matrix
presence = reshape(presence', mesh.Nx, mesh.Ny)';
O_b = reshape(O_b', mesh.Nx, mesh.Ny)';
O_d = reshape(O_d', mesh.Nx, mesh.Ny)';

% Convert matrices into "image" matrices
presence = flipud(presence);
O_b = flipud(O_b);
O_d = flipud(O_d);

end
