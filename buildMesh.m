function mesh = buildMesh(Nx, Ny, dx)
% This function simply builds a basic struct that contains the mesh
% information (for a 2D mesh), for use inside other functions
%
% Usage:
%     mesh = buildMesh(Nx, Ny, dx)
%
% INPUTS:
%
%     Nx: Number of points in the x-direction
%     Ny: Number of points in the y-direction
%     dx: Mesh spacing (assumed square)
%
% OUPUTS:
%
%   mesh: A struct with fields "points" (point locations) and "Nx" and "Ny"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a set of points that fall in centres of pixels
xv = linspace(dx/2, dx*(Nx-1/2), Nx);
yv = linspace(dx/2, dx*(Ny-1/2), Ny);
[Y,X] = meshgrid(yv,xv);
points = [X(:) Y(:)];

% Store in mesh struct
mesh.points = points;
mesh.Nx = Nx;
mesh.Ny = Ny;

end

