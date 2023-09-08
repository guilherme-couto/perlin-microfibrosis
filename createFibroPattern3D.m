function [presence, O_b, O_d, F] = createFibroPattern3D(mesh, density, params, Ps, offsets, offsets2)
% This function takes an input three-dimensional mesh of points and creates
% a three-dimensional realisation (presence/absence of deposited material)
% of fibrosis.
%
% Usage:     presence = createFibroPattern3D(mesh, density, params, Ps, offsets, offsets2)
%
% INPUTS:   mesh:      a mesh object (use createMesh.m to create)
%           density:   the proportion of points to be occupied by fibrosis
%           params:    a set of parameters defining the properties of the noise pattern (see below)
%           Ps:        an m x n matrix, which is m rows of the numbers 0:n-1 arranged in random order (permutation tables for random assignment of vectors in Perlin noise)
%           offsets:   an m x 2 matrix that specifies grid offsets for each octave in octave noise
%           offsets2:  an m x 2 matrix that specifies grid offsets for each octave in octave noise
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
%           (F):        the sinusoidal field used for fibre-type patterns (optional)
%
% PARAMS:   Parameters are provided as a single vector, specified as:
%           [ fibreness, fibre_separation_y, fibre_separation_z, patchiness, feature_size, roughness, patch_size, alignment_ratio_y, alignment_ratio_z, phi, theta ]
%
%           (phi and theta are the two angles that define the fibre direction in three dimensions)
%           (properties with an appended _y or _z refer to that property in that direction, before rotation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


%%% INITIAL SETUP

% Read out paramaters from vector
params = num2cell(params);
[fibreness, fibre_sep1, fibre_sep2, patchiness, feature_size, roughness, patch_size, fibre_alignment_y, fibre_alignment_z, phi, theta] = deal(params{:});

% Create rotated points using the 3D rotation matrix
R = [ cos(phi) * cos(theta), sin(phi), -cos(phi) * sin(theta);
      -cos(phi) * sin(theta), cos(phi), sin(phi) * sin(theta);
      sin(theta), 0, cos(theta) ];
R_points = R * mesh.points';  
  
% Create new permutation tables from the provided one by applying it to
% itself several times
for k = 1:size(Ps,1)
    Ps2(k,:) = Ps(k,Ps(k,:)+1);
    Ps3(k,:) = Ps(k,Ps2(k,:)+1);
    Ps4(k,:) = Ps(k,Ps3(k,:)+1);
end


%%% FIBRE-SELECTING FIELD

% Non-fibre cases are properly set to have zero fibreness
if isnan(fibreness)
    fibreness = 0;
    fibre_sep1 = 1;
    fibre_sep2 = 1;
end

% Standard non-tuned parameter values for fibres
n_fibres_similarity = 4;
wiggle_feature_length = 4;
phasefield_strength = 5;
       
% Create noisefields with these properties (used for phase modulation)
phasefield_points = [ R_points(1,:) / wiggle_feature_length; R_points(2,:) / (n_fibres_similarity * fibre_sep1); R_points(3,:) / (n_fibres_similarity * fibre_sep2) ];   
phasefield1 = Octave3D(phasefield_points, 4, 0.5, Ps, offsets);
phasefield2 = Octave3D(phasefield_points, 4, 0.5, Ps4, offsets2);

% Create the sinusoidal pattern using cos(RX^T [0;1]), with phase
% modulated using this phasefield
F = 0.5 + 0.5 * cos( 2*pi * (R_points(2,:) / fibre_sep1 + phasefield_strength * (phasefield1 - 0.5) ) ) .* cos( 2*pi * (R_points(3,:) / fibre_sep2 + phasefield_strength * (phasefield2 - 0.5) ) );
    
% Sharpen this field by powering it up many times (hard-coded for
% efficiency, but can be easily modified if a different 'sharpening factor'
% is desired). 
F = F .* F .* F .* F .* F .* F .* F .* F .* F .* F .* F .* F .* F .* F .* F;
    

%%% MAIN FIBROSIS PATTERNING FIELD

% Apply stretch to create stretched input points, and evaluate the
% noisefield at those points
R_points = [ R_points(1,:) / fibre_alignment_y^(1/3) / fibre_alignment_z^(1/3); R_points(2,:) * fibre_alignment_y^(2/3) / fibre_alignment_z^(1/3); R_points(3,:) / fibre_alignment_y^(1/3) * fibre_alignment_z^(2/3) ];
O_b = Octave3D( R_points / feature_size, 4, roughness, Ps2, offsets);
    
    
%%% LARGE-SCALE DENSITY VARIATION

% Use base points (non-rotated) stretched by patch size for density
% variation field
O_d = Octave3D(mesh.points' / patch_size, 4, 0.5, Ps3, offsets);
    

%%% FIBROSIS PATTERN CREATION

% Final noisefield is the weighted combination of the component fields
noise = ( 1 - fibreness + fibreness * F ) .* O_b + patchiness * O_d;

% Convert the final noisefield into a binary pattern
presence = thresholdPattern(noise, density);