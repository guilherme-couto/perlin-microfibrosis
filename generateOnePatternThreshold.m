function pattern = generateOnePatternThreshold(params, density, seed_num, mesh)
% This function takes the provided set of parameter values, and creates the
% representation of a fibrosis pattern with the desired density.
%
% INPUTS:
%
% params - the parameter values for the generator to use (1x8 vector)
% density - density of fibrosis in the pattern to be generated
% seed_num - seed for the random number generator
% (mesh) - optionally provided mesh to specify size of pattern
%
% PARAMETER INFORMATION:
%
% 1 - FIBRENESS: The extent to which patterns exhibit long, thin fibres 
%     aligned in consistent directions
%     ::: If set to NaN, a pattern without fibres will be created :::
% 2 - FIBRE SEPARATION: The average spacing between fibres (in units
%     matching those used in input mesh
% 3 - PATCHINESS: The extent of inconsistency in pattern density (high
%     patchiness will produce distinct regions of higher and lower density)
% 4 - FEATURE SIZE: The overall size of obstacle features in the mesh (in
%     units matching the mesh)
% 5 - ROUGHNESS: The roughness of feature edges (values from [0,1], may
%     cause issues for values of precisely 0 or 1)
% 6 - PATCH SIZE: The size of regions over which density varies (with 
%     extent according to PATCHINESS)
% 7 - FIBRE ALIGNMENT: The extent to which non-fibre features are aligned
%     to the fibre direction (i.e. extent of feature anisotropy)
% 8 - DIRECTION: An angle (in radians) defining the orientation of fibres
%     and/or feature anisotropy
%
% OUTPUTS:
%
% pattern - the generated pattern (2D matrix of 0s and 1s)

% Define a 'fibrosis' colormap
fibroclr = [[0.1, 0.5, 0.8]; [0.9, 0.5, 0.1]]; % Blue and orange

% Create the mesh if one wasn't provided (uses values from paper)
if nargin < 4
    mesh = buildMesh(250, 400, 1/136);
end

% Use the fibre-free generator if NaNs are present in input params
% vector, or if only non-fibre parameters provided, otherwise 
% use the standard generator
threshold = density;
[permute_table, offset_table] = generateTables(seed_num);
if any(isnan(params))
    [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params(3:8), permute_table, offset_table);
elseif length(params) == 6
    [presence, ~, ~] = createFibroPatternNoFibres(mesh, threshold, params, permute_table, offset_table);
else
    [presence, ~, ~, ~] = createFibroPattern(mesh, threshold, params, permute_table, offset_table);
end

pattern = presence;

end

function [permute_table, offset_table] = generateTables(seed)
% This function generates the permutation and offset tables for the provided seed.
%
% INPUTS:
%
% seed - the seed to be used for the random number generator
%
% OUTPUTS:
%
% permute_table - the permutation table generated
% offset_table - the offset table generated

% Set the seed for the random number generator
rng(seed);

% Assume a decent safe number like eight for the number of offsets
N_freqs = 8;

% Permutation tables for this seed
for j = 1:N_freqs
    permute_table(j,:) = int32(randperm(256) - 1);
end

% Offset table for this seed
offset_table = rand(N_freqs, 2) - 0.5;
    
end
