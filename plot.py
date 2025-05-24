import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Function to convert radians to degrees
def radians_to_degrees(radians):
    return radians * (180 / np.pi)

# Reading data from the file
file_path = f'./outputs/analysis/analysis.txt'  # Change to the path of your file
data = []

with open(file_path, 'r') as file:
    for line in file:
        line = line.strip()
        # Skip empty lines or lines containing "Distinct Seeds Used"
        if not line or "Distinct Seeds Used" in line:
            continue
        
        parts = line.split(', ')
        # Ensure that we have the expected number of parts
        if len(parts) < 6:
            print(f"Skipping line due to unexpected format: {line}")
            continue
        
        try:
            fibro_type = parts[0].split(': ')[1]
            density = float(parts[1].split(': ')[1]) * 100  # Convert density to percentage
            angle = float(parts[2].split(': ')[1])
            path_found_count = int(parts[3].split(': ')[1])
            total_executions = int(parts[4].split(': ')[1])
            percentage = float(parts[5].split(': ')[1].replace('%', ''))  # Keep percentage as is
            
            data.append({
                'Fibro Type': fibro_type,
                'Density': density,
                'Angle': angle,
                'Path Found Count': path_found_count,
                'Total Executions': total_executions,
                'Percentage': percentage
            })
        except (IndexError, ValueError) as e:
            print(f"Error processing line: {line}. Error: {e}")

# Creating a DataFrame
df = pd.DataFrame(data)

# Defining the fibro types to plot
fibro_types = df['Fibro Type'].unique()

# Creating the plots
for fibro_type in fibro_types:
    plt.figure()
    
    # Filtering data by fibro type
    subset = df[df['Fibro Type'] == fibro_type]
    
    # Creating lines for each angle
    for angle in subset['Angle'].unique():
        angle_data = subset[subset['Angle'] == angle]
        plt.plot(angle_data['Density'], angle_data['Percentage'], marker='o', label=f"{int(radians_to_degrees(angle))}°")
    
    plt.title(f"Fibro Type: {fibro_type}")
    plt.xlabel("Density (%)")
    plt.ylabel("Percentage of Path Found")
    plt.xticks(np.arange(0, 101, 10))
    plt.yticks(np.arange(0, 101, 10))
    plt.ylim(0, 100)
    plt.xlim(0, 100)
    plt.legend(title="Angles")
    plt.grid()
    
    # Saving the plot
    plt.savefig(f'{fibro_type}_graph.png', )
    plt.close()

# Creating a combined plot
plt.figure(figsize=(12, 8))

# Plotting each fibro type in the combined figure
for i, fibro_type in enumerate(fibro_types, 1):
    plt.subplot(2, 2, i)  # Create a subplot for each fibro type

    # Filtering data by fibro type
    subset = df[df['Fibro Type'] == fibro_type]
    
    for angle in subset['Angle'].unique():
        angle_data = subset[subset['Angle'] == angle]
        plt.plot(angle_data['Density'], angle_data['Percentage'], marker='o', label=f"{int(radians_to_degrees(angle))}°")
    
    plt.title(f"Fibro Type: {fibro_type}")
    plt.xlabel("Density (%)")
    plt.ylabel("Percentage of Path Found")
    plt.xticks(np.arange(0, 101, 10))
    plt.yticks(np.arange(0, 101, 10))
    plt.ylim(0, 100)
    plt.xlim(0, 100)
    plt.legend(title="Angles")
    plt.grid()

# Adjust layout and save the combined plot
plt.tight_layout()
plt.savefig('combined_fibro_types_graph.png')
plt.close()

print("Plots generated and saved successfully.")
