function problem = createProblemFromMap(occ_map, phi)
% This function creates a problem object (as used by the monodomain-solving
% and homogenisation code) from the given input map (a matrix of logicals,
% or 0's and 1's). No stimulus is set up at current, because this function
% is used together with the homogenisation code and not actually simulated.
%
% The second, optional argument specifies the angle (away from horizontal)
% that fibres should be oriented

% Define the diffusion tensor
D = [ 3, 0; 0, 1 ];    % Fibre-biased conduction (fibres horizontal)

% Define the grid separation
mesh_separation = 1/137 * 1000;   % Pixel size converted to micrometres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set rotation to zero if not requested
if nargin < 2
    phi = 0;
end

% Create a rotation matrix generator
R = @(phi) [cos(phi), sin(phi); -sin(phi), cos(phi)];

% Use rotation matrices to adjust conductivity so fibres are oriented at
% the requested angle
D = R(-phi) * D * R(phi);


% Convert the grid separation to centimetres for consistent units
mesh_separation = mesh_separation / 10000;

% Read out the size
[Ny,Nx] = size(occ_map);

% Find the domain size
Lx = Nx * mesh_separation;
Ly = Ny * mesh_separation;


%%% Define volumes of each element
Vfrac = double(~occ_map);      % Volume fraction of one in all non-occupied elements (no elements are partially occupied)


%%% Create stimulus sites

% Nodes will be placed at element boundaries (vertex-centred finite volume)
% So first create node positions
[nodeX, nodeY] = meshgrid( linspace(0,Lx,Nx+1), linspace(0,Ly,Ny+1) );

% Initialise stimulus matrix to zeroes
stim_sites1 = false(size(nodeY));
stim_sites2 = false(size(nodeY));


%%% Specify the cell model to use at all sites

% List cell models that will be used here
cell_models = {'TT3epi'};
% Assign models to cells (by number)
model_assignments = ones(size(nodeX));



%%% Process and save all data

% Read out base diffusivity levels from the diffusion tensor
D_xx = D(1,1);
D_xy = D(1,2);
D_yy = D(2,2);

% Create matrices of diffusion values, with zero in blocked regions
D_xx = D_xx * (~occ_map);
D_xy = D_xy * (~occ_map);
D_yy = D_yy * (~occ_map);

% Store problem details in the 'problem' structure
problem.occ_map = occ_map;
problem.D_tensor.D_xx = D_xx;
problem.D_tensor.D_xy = D_xy;
problem.D_tensor.D_yy = D_yy;
problem.Vfrac = Vfrac;
problem.grid.dx = mesh_separation;
problem.grid.dy = mesh_separation;
problem.grid.Lx = Lx;
problem.grid.Ly = Ly;
problem.Nx = Nx;
problem.Ny = Ny;
nodeX = nodeX'; nodeX = nodeX(:);
nodeY = nodeY'; nodeY = nodeY(:);
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.stim_sites1 = stim_sites1;
problem.stim_sites2 = stim_sites2;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

end