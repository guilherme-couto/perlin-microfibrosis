function FIGURE_showConductivity(fibro_type)
% This function generates a large slice of the requested type of
% patterning, and plots a summary of its diffusive properties by
% homogenising sub-elements of it.
%
% Fibrosis types are:
%    1 - Interstitial
%    2 - Compact
%    3 - Diffuse
%    4 - Patchy

% Specify the number of fibrosis pixels to homogenise, and their width
dx = 1/136;     % Width (in mm) of each pixel
Nh = 200;       % 15 x 1/136 = 367 micron width

% Specify the number of large elements to use in building the slice
Ex = 6;        % Number of elements in horizontal direction
Ey = 8;        % Number of elements in vertical direction

% Specify the seed number (for repeatability)
seed_num = 7;

% Specify how the size of arrows are scaled
arrow_scaling = 28;

% Fibrosis colours from the previous figures
fibro_clrs = [ [0.95, 0.85, 0.55];     % Tissue
               [0.80, 0.20, 0.20]  ];  % Fibrosis
% Darken and slightly de-saturate the old fibrosis colours
fibro_clrs = fibro_clrs * 0.6 - 0.2 * ( fibro_clrs - mean(fibro_clrs,2) );
% Colour both directions of conduction in a light blue that stands out
% against the darkened fibrosis colours
arrow_clrs = 1/255 * [ 105, 205, 255;             % Dominant direction of conduction
                       105, 205, 255 ];           % Secondary direction of conduction      

% Specify the names of fibrosis types (for use in figure titling)
fibro_types = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};
                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the generator tunings
load('param_modes.mat','param_modes','densities');

% Load in the seed information
load('fibro_seedinfo.mat','permute_tables','offset_tables');

% Generate the slice of pattern for the requested fibrosis type
gen_params = param_modes(fibro_type,:);
phi = gen_params(end);                           % Read out fibre direction
density = densities(fibro_type);

% Prepare the mesh
Nx = Ex * Nh;
Ny = Ey * Nh;
mesh = buildMesh(Nx,Ny,dx);

% Generate the pattern for the current class and seed number
if any(isnan(gen_params))
    [slice_occ, ~, ~] = createFibroPatternNoFibres(mesh, density, gen_params(~isnan(gen_params)), permute_tables{seed_num}, offset_tables{seed_num});
else
    [slice_occ, ~, ~] = createFibroPattern(mesh, density, gen_params, permute_tables{seed_num}, offset_tables{seed_num});
end     

% Convert this slice into a problem array as used by homogenisation code
problem = createProblemFromMap(flipud(slice_occ), phi);
            
% Homogenise this whole problem - use linear BCs as per previous work
homog_problem = homogeniseFull2DProblem(problem, Nh, Nh, 'linear');

% Read out homogenised problem size
[Nhy, Nhx] = size(homog_problem.D_tensor.D_xx);

% Initialise figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
hold on;

% Display the map of fibrosis but in the greyed-out colours
imagesc(problem.occ_map);
colormap(fibro_clrs);

% Plot the gridlines showing the edges of each element
for i = 1:Nhx+1
    plot([(i-1)*Nh, (i-1)*Nh],[-1 Ny+1],'k','LineWidth',2);
end
for j = 1:Nhy+1
    plot([-1 Nx+1],[(j-1)*Nh, (j-1)*Nh],'k','LineWidth',2);
end

% Display the summary arrows for each element in the homogenised problem
for i = 1:Nhx
    for j = 1:Nhy
        
        % Grab out the diffusion tensor at this point
        D = [homog_problem.D_tensor.D_xx(j,i), homog_problem.D_tensor.D_xy(j,i); homog_problem.D_tensor.D_xy(j,i), homog_problem.D_tensor.D_yy(j,i)];
        % Find its principal directions (and strengths) of diffusion
        [V, LAM] = eig(D);
        % Get the order of eigenvalues in decreasing strength
        [~,I] = sort(diag(LAM),'descend');
        
        % Find the centre of this "element"
        ele_centre = [ (2*i - 1) * Nh/2, (2*j - 1) * Nh/2 ];
        
        % Plot the two arrows
        for k = 1:2
            
            % Read out the arrow vector here
            dp = LAM(I(k),I(k)) * V(:,I(k))';
            % Flip it to point left to right if it doesn't already
            if dp(1) < 0
                dp = -dp;
            end
            % Plot the arrow
            quiver(ele_centre(1), ele_centre(2), dp(1)*arrow_scaling, dp(2)*arrow_scaling,'LineWidth',3,'Color',arrow_clrs(k,:),'MaxHeadSize',15);
        end
        
    end
end

% Ensure each pixel is actually square
axis equal;
axis off;

% Trim plot to just the area of interest
xlim([0 Nx]);
ylim([0 Ny]);

% Add title
title(fibro_types{fibro_type},'FontSize',24);