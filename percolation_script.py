from functions import *
from collections import OrderedDict

angles_degree = np.arange(0, 91, 15)
angles = [f"{np.radians(a):.4f}" for a in angles_degree] # radians
print(f"Angles: {angles}")

import random
# seed_values = random.sample(range(1, 500), 20)
seed_values = [214, 483, 454, 163, 118, 139, 182, 22, 480, 451,
               449, 141, 364, 251, 487, 87, 164, 20, 80, 21,
               496, 369, 179, 176, 281, 111, 491, 333, 7, 309,
               407, 254, 340, 314, 190, 486, 125, 18, 452, 14,
               294, 109, 47, 404, 315, 401, 326, 239, 114, 51,
               238, 395, 160, 418, 434, 84, 378, 336, 285, 92,
               11, 244, 421, 497, 278, 13, 8, 405, 50, 191,
               137, 66, 97, 471, 478, 4, 53, 183, 341, 284,
               391, 127, 217, 55, 94, 415, 236, 437, 424, 392,
               489, 36, 174, 158, 290, 30, 300, 302, 126, 275,
               37, 181, 154, 351, 145, 468, 155, 201, 54, 467]

save_plots = False
save_fibrosis_map = False
mode = 'composition'
maps_dir = f'./using_{mode}/outputs/maps_figures'
result_path = f'./using_{mode}/outputs/results/result.csv'
analysis_path = f'./using_{mode}/outputs/analysis/analysis.txt'
tex_tables_path = f'./using_{mode}/outputs/analysis/tables.tex'

# Iterate over values
# for angle in angles:
#     fibrosis_params = [ ['interstitial', f'[0.3, 0.31, 0.32, 0.24, 0.96, 4.67, 1.89, {angle}]'],
#                         ['compact',      f'[NaN, NaN, 0.44, 0.96, 0.59, 2.03, 2.47, {angle}]'],
#                         ['diffuse',      f'[NaN, NaN, 0.49, 0.07, 0.44, 1.22, 2.17, {angle}]'],
#                         ['patchy',       f'[0.38, 0.31, 0.42, 0.32, 0.78, 2.1, 2.5, {angle}]']]

#     for elem in fibrosis_params:
#         fibro_typename = elem[0]
#         params = elem[1]

#         densities = ['0.100', '0.200', '0.300', '0.400', '0.500', '0.600', '0.700', '0.800', '0.900']
        
#         if fibro_typename == 'compact':
#             densities_extension = np.linspace(0.31, 0.69, 20)
#         elif fibro_typename == 'diffuse':
#             densities_extension = np.linspace(0.41, 0.59, 10)
#         elif fibro_typename == 'interstitial':
#             densities_extension = np.linspace(0.31, 0.59, 15)
#         elif fibro_typename == 'patchy':
#             densities_extension = np.linspace(0.21, 0.59, 20)
#         densities_extension = [f"{density:.3f}" for density in densities_extension]
#         densities = list(OrderedDict.fromkeys(densities + densities_extension))

#         print(densities)
        
#         for density in densities:
#             for seed in seed_values:
#                 grid_name = f'{fibro_typename}_{density}_{seed}_{angle}'

#                 # Check if simulation is already done
#                 if is_simulation_done(result_path, fibro_typename, density, seed, angle):
#                     print(f'Simulation {grid_name} already done. Skipping...')
#                     continue

#                 print(f'\n#########################################################')
#                 print(f'PERCOLATION FOR {grid_name}\n')

#                 # Generate the mesh
#                 execute_fibrosis_generator(seed, fibro_typename, params, density)

#                 # Read the mesh (.alg) file
#                 mesh_file_path = f'./fibrosis_generator/original_mesh_alg/{fibro_typename}_density_{density}_seed_{seed}.alg'
#                 data = read_alg_file(mesh_file_path)

#                 # Create a 2D grid of elements from the data
#                 grid = create_2d_grid(data)

#                 # Create a fibrosis map from the grid
#                 fibrosis_map = create_fibrosis_map(grid)

#                 # Plot and save the fibrosis map
#                 if save_fibrosis_map or save_plots:
#                     if not os.path.exists(maps_dir):
#                         os.makedirs(maps_dir)
#                     fibrosis_map_path = os.path.join(maps_dir, f'{grid_name}.png')
#                     plot_fibrosis_map(fibrosis_map, fibrosis_map_path, grid_name)

#                 # Call percolation process
#                 path_found = percolation(grid, grid_name, save_plots=save_plots)
#                 os.remove(mesh_file_path)

#                 # Save result to csv
#                 append_result_to_csv(result_path, fibro_typename, density, seed, angle, path_found)

#                 if save_plots:
#                     # Create GIF from plots of percolation process
#                     create_gif_from_plots(grid_name)

# Process CSV file and write the analysis to a TXT and a LaTeX files
process_fibro_data(result_path, analysis_path, tex_tables_path)