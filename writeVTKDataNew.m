function writeVTKDataNew(mesh3D, presence, filename)
% This function takes the input 3D noise field, and writes a VTK file which
% creates a square surface polygon on the outside edges of each voxel
% clump of fibrosis.

% Assign default filename and file title if none given
if nargin < 3
    filename = 'presence3D';
end

% Convert the data into a 3D array
Nx = mesh3D.N_points(1);
Ny = mesh3D.N_points(2);
Nz = mesh3D.N_points(3);
P = reshape(presence, Ny,Nx,Nz);

% Pad the outside of P with zeroes
Pp = zeros(Ny+2,Nx+2,Nz+2);
Pp(2:Ny+1,2:Nx+1,2:Nz+1) = P;

% Find the grid spacings using the point data
xv = mesh3D.points(:,1);
uxv = unique(xv);
dx = uxv(2) - uxv(1);
yv = mesh3D.points(:,2);
uyv = unique(yv);
dy = uyv(2) - uyv(1);
zv = mesh3D.points(:,3);
uzv = unique(zv);
dz = uzv(2) - uzv(1);


% Create a grid of points - simply arrange the numbers 1 to N_points
% into a 3D array so that polygon vertices can be specified in terms of
% node IDs - note that node IDs start at zero
point_IDs = 0:(Nx+1)*(Ny+1)*(Nz+1)-1;
point_IDs = reshape( point_IDs, Ny+1, Nx+1, Nz+1);

[X,Y,Z] = meshgrid( min(uxv)-dx/2:dx:max(uxv)+dx/2+sqrt(eps), min(uyv)-dy/2:dy:max(uyv)+dy/2+sqrt(eps), min(uzv)-dz/2:dz:max(uzv)+dz/2+sqrt(eps) );
vertex_points = [X(:), Y(:), Z(:)];



% Find mismatches in presence between neighbouring voxels, indicating a
% surface to draw - separate lists are formed for whether it is the
% "former" or "latter" element that is occupied, so that normals can be
% correctly specified later
xfill_lopen = ~Pp(2:Ny+1,1:Nx+1,2:Nz+1) & Pp(2:Ny+1, 2:Nx+2, 2:Nz+1 );
xfill_ropen = Pp(2:Ny+1,1:Nx+1,2:Nz+1) & ~Pp(2:Ny+1, 2:Nx+2, 2:Nz+1 );
yfill_dopen = ~Pp(1:Ny+1,2:Nx+1,2:Nz+1) & Pp(2:Ny+2, 2:Nx+1, 2:Nz+1 );
yfill_uopen = Pp(1:Ny+1,2:Nx+1,2:Nz+1) & ~Pp(2:Ny+2, 2:Nx+1, 2:Nz+1 );
zfill_iopen = ~Pp(2:Ny+1,2:Nx+1,1:Nz+1) & Pp(2:Ny+1, 2:Nx+1, 2:Nz+2 );
zfill_oopen = Pp(2:Ny+1,2:Nx+1,1:Nz+1) & ~Pp(2:Ny+1, 2:Nx+1, 2:Nz+2 );

% Create a list of sites to fill (a cell array because each will have its
% own properties
fill_sites = { xfill_lopen, xfill_ropen, yfill_dopen, yfill_uopen, zfill_iopen, zfill_oopen };

% Specify the indices used to map to co-ordinates of the voxels where
% filling is to occur
coord_inds = { [Ny, Nx+1, Nz], [Ny, Nx+1, Nz], [Ny+1, Nx, Nz], [Ny+1, Nx, Nz], [Ny, Nx, Nz+1], [Ny, Nx, Nz+1] };

% Specify index shifts for the three face corners away from the base corner
% (simply specifying how index shifts as we step in x, y, or z direction as
% appropriate)
x_shifts = [1, (Ny+1)*(Nx+1)+1, (Ny+1)*(Nx+1),];
y_shifts = [Ny+1, (Ny+1)*(Nx+1)+(Ny+1), (Ny+1)*(Nx+1)];
z_shifts = [1, Ny+1+1, Ny+1];
ind_shifts = { x_shifts, x_shifts, y_shifts, y_shifts, z_shifts, z_shifts };

% Specify the outward-facing normals for faces in each direction ( e.g.
% left-open face corresponds to a normal of [-1, 0, 0] )
dir_normals = { [-1 0 0], [1 0 0], [0 -1 0], [0 1 0], [0 0 -1], [0 0 1] };

% Loop over the whole set of lists of faces to fill in
polygon_vertices = [];
polygon_normals = [];
for k = 1:length(fill_sites)
    
    % Read out the sites
    foundsites = find( fill_sites{k} );
    % Convert to co-ordinates
    [fill_i, fill_j, fill_k] = ind2sub( coord_inds{k}, foundsites );
    % Convert co-ordinates back to locations in the grid of points
    % (this is different to the locations in foundsites)
    fill_ind = sub2ind([Ny+1, Nx+1, Nz+1], fill_i, fill_j, fill_k);
    
    % Append the vertices of this list of faces to the face vertices list
    polygon_vertices = [polygon_vertices; [ point_IDs( fill_ind ), point_IDs( fill_ind + ind_shifts{k}(1) ), point_IDs( fill_ind + ind_shifts{k}(2) ), point_IDs( fill_ind + ind_shifts{k}(3) ) ] ];

    % Append the normals information for these faces
    polygon_normals = [polygon_normals; repmat( dir_normals{k}, length(fill_ind), 1 ) ];
    
end

%%% Remove any nodes that aren't involved in polygons

% Find list of vertices to keep
[poly_vertexlist,~,IC] = unique(polygon_vertices(:)+1);
% Grab out only the kept vertices as the new vertices
vertex_points = vertex_points(poly_vertexlist,:);
Npts = size(vertex_points,1);
% Create a mapping onto [0,Npts-1]
mapping = NaN( (Nx+1)*(Ny+1)*(Nz+1), 1 );
mapping(poly_vertexlist) = (1:length(poly_vertexlist)) - 1;
% Relabel all vertices in the polygon list with their new IDs using the mapping
polygon_vertices = reshape( mapping(poly_vertexlist(IC)), size(polygon_vertices) );


% Complete the polygon vertex write data by appending the number of
% vertices as the first column
N_polys = size(polygon_vertices,1);
polygon_writedata = [ 4*ones( N_polys, 1), polygon_vertices ];


% Find centre of mass of each solid element
Ex = reshape( vertex_points(polygon_vertices+1,1), size(polygon_vertices) );
Ey = reshape( vertex_points(polygon_vertices+1,2), size(polygon_vertices) );
Ez = reshape( vertex_points(polygon_vertices+1,3), size(polygon_vertices) );

CoM = [ sum(Ex,2)/4, sum(Ey,2)/4, sum(Ez,2)/4 ];

% Open the text file
fobj = fopen([filename,'.vtk'],'wt');

% Write header lines
fprintf(fobj, '# vtk DataFile Version 2.0\n');
fprintf(fobj, '%s\n', filename);
fprintf(fobj, 'ASCII\n');
fprintf(fobj, 'DATASET POLYDATA\n');
fprintf(fobj, 'POINTS %u float\n', Npts);
fprintf(fobj, '%g %g %g\n', vertex_points');

% Write the polygon data
fprintf(fobj, 'POLYGONS %u %u\n', N_polys, N_polys*5);
fprintf(fobj, '%g %g %g %g %g\n', polygon_writedata');

% Write the normals data
fprintf(fobj, 'CELL_DATA %g\n', N_polys);
fprintf(fobj, 'NORMALS cell_normals float\n');
fprintf(fobj, '%g %g %g\n', polygon_normals');

% Write the scalar data used for colouration
fprintf(fobj, 'SCALARS cell_scalars float 1\n');
fprintf(fobj, 'LOOKUP_TABLE default\n');
fprintf(fobj, '%g\n', CoM(:,3));

% Close the text file
fclose(fobj);

end

