% Script to generate fibrosis patterns and compute ellipse-based metrics
% for four types of fibrosis: compact, diffuse, interstitial, and patchy.
% Each type generates 100 patterns using specified parameters and density.
% Metrics are saved into CSV files in the specified format.

clear; clc;

% Add any required paths here (e.g., for functions)
% addpath('path/to/generator');

% Define a 'fibrosis' colormap
% fibroclr = [[0.95, 0.85, 0.55]; [0.8, 0.2, 0.2]]; % Brown and red
fibroclr = [[0.1, 0.5, 0.8]; [0.9, 0.5, 0.1]];

% Power thresholds for ellipseMetrics
power_thresholds = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

% Fibrosis types and parameters
fibrosis_types = {'interstitial', 'compact', 'diffuse', 'patchy'};
fibrosis_params = {
    [0.3, 0.31, 0.32, 0.24, 0.96, 4.67, 1.89, 0.59341], 0.096;
    [NaN, NaN, 0.44, 0.96, 0.59, 2.03, 2.47, -0.15708], 0.472;
    [NaN, NaN, 0.49, 0.07, 0.44, 1.22, 2.17, 0.19199], 0.22;
    [0.38, 0.31, 0.42, 0.32, 0.78, 2.1, 2.5, 1.18680], 0.269
};

% Common tolerance for density matching
tolerance = 0.005;

% Loop over each fibrosis type
for t = 1:numel(fibrosis_types)
    type_name = fibrosis_types{t};
    params = fibrosis_params{t,1};
    density = fibrosis_params{t,2};

    % Generate header names in the order: orientation_X, major_axis_X, minor_axis_X
    thresholds = 10:10:90;
    header_names = {'fibro_typename', 'seed'};
    for i = 1:length(thresholds)
        th_str = num2str(thresholds(i));
        header_names = [header_names, ...
            strcat('orientation_', th_str), ...    
            strcat('major_axis_', th_str), ...
            strcat('minor_axis_', th_str)];
    end

    % Initialize results cell array (1 header row + 100 data rows, each with 29 columns - fibro_typename, seed, 27 metrics))
    results = cell(101, length(header_names));
    results(1, :) = header_names;

    % Collect 9 patterns for plotting
    example_patterns = {};

    for i = 1:100
        fprintf('Generating %s pattern %d/100\n', type_name, i);
        seed = i; % could be randomized
        pattern = generateOnePatternComposition(params, density, seed, tolerance);

        % Save first 9 for visualization
        if i <= 9
            example_patterns{end+1} = pattern;
        end

        % Compute metrics
        metrics = ellipseMetrics(pattern, power_thresholds);

        % Store results
        results{i+1,1} = type_name;
        results{i+1,2} = seed;
        results(i+1,3:end) = num2cell(metrics);
    end

    % Save CSV
    filename = sprintf('metrics_%s.csv', type_name);
    fid = fopen(filename, 'w');

    % Write header
    for j = 1:numel(results(1,:))-1
        fprintf(fid, '%s,', results{1,j});
    end
    fprintf(fid, '%s\n', results{1,end});

    % Write data rows
    for i = 2:size(results,1)
        fprintf(fid, '%s,%d', results{i,1}, results{i,2});
        for j = 3:size(results,2)
            fprintf(fid, ',%.6f', results{i,j});
        end
        fprintf(fid, '\n');
    end
    fclose(fid);

    % Plot 9 example patterns (3x3 grid)
    figure('visible','off');
    colormap(fibroclr);
    for k = 1:9
        subplot(3,3,k);
        imagesc(example_patterns{k});
        axis equal off;
        title(sprintf('%s #%d', type_name, k));
    end
    img_name = sprintf('examples_%s.png', type_name);
    print(img_name, '-dpng');
    close;
end