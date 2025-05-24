import numpy as np
import matplotlib.pyplot as plt
import os, re, csv, subprocess
import pandas as pd

class Element:
    def __init__(self, fibrosis_tag, activated=False):
        self.fibrosis = fibrosis_tag == 1
        self.activated = activated

def read_alg_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    data = [list(map(float, line.strip().split(','))) for line in lines]
    
    print(f'File {file_path} read successfully. Number of elements: {len(data)*len(data[0])}')
    return data

def create_2d_grid(data):
    # Get unique x and y values
    x_values = sorted(set(row[0] for row in data))
    y_values = sorted(set(row[1] for row in data))

    # Create an empty 2D grid
    grid = np.empty((len(y_values), len(x_values)), dtype=object)

    # Map the data to the grid using the x and y values as indices and fibrrosis tag as the value
    for row in data:
        center_x, center_y, _, _, _, _, fibrosis_tag, _, _, _ = row
        i = y_values.index(center_y)
        j = x_values.index(center_x)
        grid[i][j] = Element(int(fibrosis_tag))

    print(f'Grid created with shape {grid.shape}')
    return grid

def grid_reached_last_column(grid):
    # Check if any element in the last column was activated
    last_column_index = len(grid[0]) - 1
    for i in range(len(grid)):
        if grid[i][last_column_index].activated:
            return True
    return False

def create_fibrosis_map(grid):
    fibrosis_map = np.zeros((len(grid), len(grid[0])), dtype=float)

    for i in range(len(grid)):
        for j in range(len(grid[i])):
            element = grid[i][j]
            if element:
                if element.fibrosis:
                    fibrosis_map[i][j] = 1  # Fibrosis
                elif element.activated:
                    fibrosis_map[i][j] = 0.5  # Activated
                else:
                    fibrosis_map[i][j] = 0  # Non fibrosis and non activated

    return fibrosis_map

def plot_fibrosis_map(fibrosis_map, output_file, title=''):
    plt.figure(figsize=(5, 5))
    plt.imshow(fibrosis_map, cmap='viridis', origin='lower')
    plt.axis('off')
    if title != '':
        plt.title(title)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    if title != '':
        print(f'Fibrosis map saved to {output_file}')

def natural_sort_key(filename):
    # Key for natural sorting of filenames
    return [int(s) if s.isdigit() else s for s in re.split('(\d+)', filename)]

def create_gif_from_plots(grid_name):
    # Check if plots_dir exists
    plots_dir = f'./outputs/gifs/percolation_plots/{grid_name}'
    if not os.path.exists(plots_dir):
        print(f"Error: Directory '{plots_dir}' not found, GIF wont be created!")
        return

    # Create a GIF from the plots
    import imageio.v2 as imageio
    save_dir = f'./outputs/gifs'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    gif_path = os.path.join(save_dir, f'percolation_{grid_name}.gif')

    # To save all the images
    images = []
    
    # Read and sort the filenames in the directory
    frames_count = 0
    for filename in sorted(os.listdir(plots_dir), key=natural_sort_key):
        if filename.endswith('.png'):
            file_path = os.path.join(plots_dir, filename)
            images.append(imageio.imread(file_path))
            frames_count += 1
    imageio.mimsave(gif_path, images, duration=0.1)
    print(f'GIF created with {frames_count} frames and saved to {gif_path}')

    # Clean directory
    for filename in os.listdir(plots_dir):
        if filename.endswith('.png'):
            file_path = os.path.join(plots_dir, filename)
            os.remove(file_path)    

def execute_fibrosis_generator(seed_num, fibro_typename, params, density):
    command = f"octave --no-gui --no-window-system --silent ./script_with_args.m -- {seed_num} '{fibro_typename}' '{params}' {density}"
    print(f'Running command {command}')
    subprocess.run(command, shell=True, cwd='fibrosis_generator')

def percolation(grid, grid_name, save_plots=True):
    # Percolation process
    rows, cols = len(grid), len(grid[0])
    iter_count = 0

    if save_plots:
        # Create a directory to save the plots
        plots_dir = f'./outputs/gifs/percolation_plots/{grid_name}'
        if not os.path.exists(plots_dir):
            os.makedirs(plots_dir)

    # List to keep track of the activated elements in the current iteration
    current_activated = []

    # Initial activation: activate the elements in the first column that are not fibrosis
    for i in range(rows):
        if not grid[i][0].fibrosis:
            grid[i][0].activated = True
            current_activated.append((i, 0))  # Add to the list of activated

    if save_plots:
        # Plot the initial state
        fibrosis_map = create_fibrosis_map(grid)
        plot_path = os.path.join(plots_dir, f'iteration_{iter_count}.png')
        plot_fibrosis_map(fibrosis_map, plot_path)

    while not grid_reached_last_column(grid):
        new_activation = False
        next_to_activate = []  # List of elements to activate in the next iteration

        # Iterate over the activated elements in the previous iteration
        for i, j in current_activated:
            neighbors = []

            if i > 0:  # Upper
                neighbors.append((i-1, j))
            if i < rows - 1:  # Bottom
                neighbors.append((i+1, j))
            if j > 0:  # Left
                neighbors.append((i, j-1))
            if j < cols - 1:  # Right
                neighbors.append((i, j+1))

            # Activate the non-activated neighbors that are not fibrosis
            for ni, nj in neighbors:
                if not grid[ni][nj].fibrosis and not grid[ni][nj].activated:
                    grid[ni][nj].activated = True
                    next_to_activate.append((ni, nj))  # Add to the list for the next iteration
                    new_activation = True
        
        # Update the list of activated elements
        current_activated = next_to_activate

        # Increment the iteration count
        iter_count += 1
        
        if save_plots:
            # Plot every 4 iterations
            if iter_count % 4 == 0:
                fibrosis_map = create_fibrosis_map(grid)
                plot_path = os.path.join(plots_dir, f'iteration_{iter_count}.png')
                plot_fibrosis_map(fibrosis_map, plot_path)

        # Break if no new activation was made
        if not new_activation:
            print(f'At iter {iter_count}, no new activation was made.')
            break

    # Create a directory to save the last plot and save it
    save_last_dir = './outputs/last_plots'
    if not os.path.exists(save_last_dir):
        os.makedirs(save_last_dir)
    fibrosis_map = create_fibrosis_map(grid)
    plot_last_path = os.path.join(save_last_dir, f'last_plot_{grid_name}.png')
    # plot_fibrosis_map(fibrosis_map, plot_last_path, title=f'At iteration {iter_count}')

    if save_plots:
        # Save the last plot (for gif)
        fibrosis_map = create_fibrosis_map(grid)
        plot_path = os.path.join(plots_dir, f'iteration_{iter_count}.png')
        plot_fibrosis_map(fibrosis_map, plot_path)

    # Verify if the grid reached the last column
    if grid_reached_last_column(grid):
        print(f'Path found from left to right for {grid_name} in {iter_count} iterations.')
        return True
    print(f'Path not found for {grid_name}.')
    return False

def append_result_to_csv(csv_path, fibro_typename, density, seed, angle, path_found):
    # Ensure the directory of the csv_path exists
    directory = os.path.dirname(csv_path)
    if not os.path.exists(directory):
        print(f"Notice: The directory {directory} does not exist. Creating it now.")
        os.makedirs(directory)  # Create the directory structure

    # Check if the file path exists
    if not os.path.exists(csv_path):
        print(f"Notice: The file {csv_path} does not exist. Creating it now.")

        # Create the file and write the header
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['fibro_typename', 'density', 'seed', 'angle', 'path_found'])
        print("File created and header written.")

    # Define the expected header
    header = ['fibro_typename', 'density', 'seed', 'angle', 'path_found']
    
    # Check if the file already exists and if it has the correct header
    file_exists = os.path.isfile(csv_path)
    if file_exists:
        with open(csv_path, 'r') as f:
            existing_header = f.readline().strip().split(',')
            if existing_header != header:
                print(f"Error: The header of the file does not match the expected. Nothing will be written. Expected: {header}")
                return
    else:
        # If the file does not exist, create it and write the header
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)

    # Append the information to the file
    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([fibro_typename, density, seed, angle, 1 if path_found else 0])
        print(f"Result appended to {csv_path} successfully.")


def is_simulation_done(result_path, fibro_typename, density, seed, angle):
    try:
        with open(result_path, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if (row['fibro_typename'] == fibro_typename and
                    row['density'] == density and
                    row['seed'] == str(seed) and
                    row['angle'] == angle):
                    return True
    except FileNotFoundError:
        # If file does not exist, return False
        pass
    return False

def process_fibro_data(input_csv, output_txt, output_latex):
    # Read the CSV file
    df = pd.read_csv(input_csv)

    # Create a DataFrame with distinct seeds
    distinct_df = df.drop_duplicates(subset=['fibro_typename', 'density', 'seed', 'angle'])

    # Group by fibro_typename, density, and angle
    grouped = distinct_df.groupby(['fibro_typename', 'density', 'angle'])

    results = []
    all_seeds = distinct_df['seed'].unique()  # Collect all distinct seeds

    previous_fibro_typename = None
    previous_density = None

    for (fibro_typename, density, angle), group in grouped:
        # Count executions where path_found is 1
        path_found_count = group['path_found'].sum()
        total_executions = len(group)

        # Calculate percentage
        percentage = (path_found_count / total_executions) * 100 if total_executions > 0 else 0

        # Store the result
        results.append({
            'fibro_typename': fibro_typename,
            'density': density,
            'angle': angle,
            'path_found_count': path_found_count,
            'total_executions': total_executions,
            'percentage': percentage
        })

    # Ensure the directory of the output_txt exists
    directory_txt = os.path.dirname(output_txt)
    if not os.path.exists(directory_txt):
        print(f"Notice: The directory {directory_txt} does not exist. Creating it now.")
        os.makedirs(directory_txt)
    
    # Ensure the directory of the output_latex exists
    directory_latex = os.path.dirname(output_latex)
    if not os.path.exists(directory_latex):
        print(f"Notice: The directory {directory_latex} does not exist. Creating it now.")
        os.makedirs(directory_latex)

    # Write the results to a TXT file
    with open(output_txt, 'w') as f:
        for result in results:
            if previous_fibro_typename != result['fibro_typename']:
                if previous_fibro_typename is not None:  # Don't add a blank line before the first entry
                    f.write('\n\n\n')
                previous_fibro_typename = result['fibro_typename']

            if previous_density != result['density']:
                if previous_density is not None:
                    f.write('\n')
                previous_density = result['density']

            f.write(f"Fibro Type: {result['fibro_typename']}, Density: {result['density']}, "
                     f"Angle: {result['angle']}, Path Found Count: {result['path_found_count']}, "
                     f"Total Executions: {result['total_executions']}, "
                     f"Percentage: {result['percentage']:.2f}%\n")
    
        # Write distinct seeds at the end of the file
        f.write(f'\nDistinct Seeds Used (total {len(all_seeds)}): ')
        for seed in all_seeds:
            f.write(f"{seed}, ")
        f.write('\n')
    
    print(f"Analysis written to {output_txt} successfully.")

    # Create a structure to group results by fibro_typename and angle
    grouped_results = {}
    for result in results:
        key = (result['fibro_typename'], result['angle'])
        if key not in grouped_results:
            grouped_results[key] = []
        grouped_results[key].append(result)

    # Write the results to a LaTeX file
    with open(output_latex, 'w') as f:
        for (fibro_typename, angle), group in grouped_results.items():
            angle_degree = int(np.degrees(float(angle)))
            f.write(rf"For \textit{{{fibro_typename}}} with an angle of {angle} radians ({angle_degree}°), the results are shown in Table \ref{{tab:percolation-{fibro_typename}-{angle_degree}}}." + '\n')
            f.write(r"\begin{table}[h!]" + '\n')
            f.write(r"    \centering" + '\n')
            f.write(rf"    \caption{{Percolation Study for \textit{{{fibro_typename}}} with Angle {angle} radians ({angle_degree}°) and 110 Different Seeds.}}" + '\n')
            f.write(r"    \begin{tabular}{|c|c|c|}" + '\n')
            f.write(r"        \hline" + '\n')
            f.write(r"        \textbf{Density} & \textbf{Paths Found} & \textbf{Percentage} \\" + '\n')
            f.write(r"        \hline" + '\n')

            # Write each result for the current fibro_typename and angle
            for r in group:
                f.write(f"        {r['density']} & {r['path_found_count']} & {r['percentage']:.2f}\\% \\\\ \hline" + '\n')

            f.write(r"    \end{tabular}" + '\n')
            f.write(rf"    \label{{tab:percolation-{fibro_typename}-{angle_degree}}}" + '\n')
            f.write(r"\end{table}" + '\n\n')
    
    print(f"Tables written to {output_latex} successfully.")

    # Plot data
    plot_data(output_txt)

def plot_data(txt_file_path):
    # Check if file exists
    if not os.path.exists(txt_file_path):
        print(f"Error: Directory '{txt_file_path}' not found, plots wont be created!")
        return

    # Reading data from the file
    data = []
    with open(txt_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Skip empty lines or lines containing "Distinct Seeds Used"
            if not line or "Distinct Seeds Used" in line:
                continue
            
            parts = line.split(', ')
            # Ensure that we have the expected number of parts
            if len(parts) != 6:
                print(f"Skipping line due to unexpected format: {line}")
                continue
            
            fibro_type = parts[0].split(': ')[1]
            density = float(parts[1].split(': ')[1]) * 100  # Convert density to percentage
            angle = float(parts[2].split(': ')[1])
            path_found_count = int(parts[3].split(': ')[1])
            total_executions = int(parts[4].split(': ')[1])
            percentage = float(parts[5].split(': ')[1].replace('%', ''))  # Keep percentage as is
            
            # Only angles greater or equal to 0.0 radians
            if angle < 0.0:
                continue
            
            data.append({
                'Fibro Type': fibro_type,
                'Density': density,
                'Angle': angle,
                'Path Found Count': path_found_count,
                'Total Executions': total_executions,
                'Percentage': percentage
            })

    # Creating a DataFrame
    df = pd.DataFrame(data)

    # Defining the fibro types to plot
    fibro_types = df['Fibro Type'].unique()

    # Ensure the directory of the output_txt exists
    dir = 'percolation_plots/'
    if not os.path.exists(dir):
        print(f"Notice: The directory {dir} does not exist. Creating it now.")
        os.makedirs(dir)

    # Creating the plots
    for fibro_type in fibro_types:
        plt.figure(figsize=(6, 5))
        
        # Filtering data by fibro type
        subset = df[df['Fibro Type'] == fibro_type]
        
        # Creating lines for each angle
        for angle in subset['Angle'].unique():
            angle_data = subset[subset['Angle'] == angle]
            plt.plot(angle_data['Density'], angle_data['Percentage'], marker='.', label=f"{int(np.degrees(float(angle)))}°")
        
        plt.title(f"{fibro_type.capitalize()}", fontweight='bold')
        plt.xlabel("Fibrosis Density (%)")
        plt.ylabel("% Path Found")
        plt.xticks(np.arange(0, 101, 10))
        plt.yticks(np.arange(0, 101, 10))
        plt.ylim(0, 100)
        plt.xlim(0, 100)
        plt.legend(title="Angles")
        plt.grid()
        
        # Saving the plot
        plt.savefig(f'percolation_plots/{fibro_type}_graph.svg')
        plt.savefig(f'percolation_plots/{fibro_type}_graph.png')
        plt.close()

    # Creating a combined plot
    # plt.figure(figsize=(12, 10))
    plt.figure(figsize=(32, 8))

    # Plotting each fibro type in the combined figure
    for i, fibro_type in enumerate(fibro_types, 1):
        plt.subplot(1, 4, i)  # Create a subplot for each fibro type

        # Filtering data by fibro type
        subset = df[df['Fibro Type'] == fibro_type]
        
        # Creating lines for each angle
        for angle in subset['Angle'].unique():
            angle_data = subset[subset['Angle'] == angle]
            plt.plot(angle_data['Density'], angle_data['Percentage'], label=f"{int(np.degrees(float(angle)))}°")
        
        plt.title(f"{fibro_type.capitalize()}", fontweight='bold', fontsize=22)
        plt.xlabel("Fibrosis Density (%)", fontsize=20)
        plt.ylabel("% Path Found", fontsize=20)
        plt.xticks(np.arange(0, 101, 10))
        plt.yticks(np.arange(0, 101, 10))
        plt.ylim(0, 100)
        plt.xlim(0, 100)
        plt.legend(fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=18)
        plt.grid()

    # Adjust layout and save the combined plot
    plt.tight_layout()
    plt.savefig('percolation_plots/combined_fibro_types_graph_horizontal.svg')
    plt.savefig('percolation_plots/combined_fibro_types_graph_horizontal.png')
    plt.close()

    print("Plots generated and saved successfully in percolation_plots.")
