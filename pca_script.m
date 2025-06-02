% Script to generate and analyze fibrosis patterns using ellipse-based metrics.
% Supports two generation modes: 'threshold' and 'composition'.
% Iterates over multiple theta values and densities for each fibrosis type.
% Saves all metrics in a single CSV file and saves 9 sample images per mode/type combo.

clear; clc;

% Define colormap
fibroclr = [[0.1, 0.5, 0.8]; [0.9, 0.5, 0.1]];

% Power thresholds for ellipseMetrics
power_thresholds = 0.1:0.1:0.9;
threshold_labels = 10:10:90;

% Output CSV initialization
csv_filename = 'fibrosis_metrics.csv';
fid = fopen(csv_filename, 'w');
fprintf(fid, 'generation_mode,fibro_typename,density,theta,seed');
for i = 1:length(threshold_labels)
    fprintf(fid, ',orientation_%d,major_axis_%d,minor_axis_%d', ...
        threshold_labels(i), threshold_labels(i), threshold_labels(i));
end
fprintf(fid, '\n');

% Fiber orientation angles
thetas_radians = [-pi/2, -pi/3, -pi/6, 0, pi/6, pi/3, pi/2];

% Fibrosis types and default parameters
fibrosis_types = {'interstitial', 'compact', 'diffuse', 'patchy'};
fibrosis_params = {
    [0.3, 0.31, 0.32, 0.24, 0.96, 4.67, 1.89, 0.59341], [0.096], [0.59341, thetas_radians];
    [NaN, NaN, 0.44, 0.96, 0.59, 2.03, 2.47, -0.15708], [0.472], [0.15708, thetas_radians];
    [NaN, NaN, 0.49, 0.07, 0.44, 1.22, 2.17, 0.19199], [0.22], [0.19199, thetas_radians];
    [0.38, 0.31, 0.42, 0.32, 0.78, 2.1, 2.5, 1.18680], [0.269], [1.18680, thetas_radians]
};

% Parameters
samples_per_config = 100;
generation_modes = {'threshold', 'composition'};
tolerance = 0.005;

% Loop over each generation mode
for g = 1:length(generation_modes)
    mode = generation_modes{g};

    % Loop over each fibrosis type
    for t = 1:numel(fibrosis_types)
        type_name = fibrosis_types{t};
        base_params = fibrosis_params{t,1};
        density_list = fibrosis_params{t,2};
        theta_list = unique(fibrosis_params{t,3});

        % Collect 9 patterns for plotting
        example_patterns = {};

        % Loop over theta and density variations
        for th = 1:length(theta_list)
            theta = theta_list(th);
            for d = 1:length(density_list)
                density = density_list(d);
                params = base_params;
                params(8) = theta; % orientation angle

                fprintf('Generating [%s] patterns for %s at theta=%.6f, density=%.3f\n', ...
                    mode, type_name, theta, density);

                for i = 1:samples_per_config
                    seed = i;

                    % Generate pattern
                    if strcmp(mode, 'threshold')
                        pattern = generateOnePatternThreshold(params, density, seed);
                    elseif strcmp(mode, 'composition')
                        pattern = generateOnePatternComposition(params, density, seed, tolerance);
                    else
                        error('Unsupported generation mode.');
                    end

                    % Save up to 9 patterns for image
                    if i <= 9 && th == 1 && d == 1
                        example_patterns{end+1} = pattern;
                    end

                    % Compute metrics
                    metrics = ellipseMetrics(pattern, power_thresholds);

                    % Write row to CSV
                    fprintf(fid, '%s,%s,%.5f,%.5f,%d', ...
                        mode, type_name, density, theta, seed);
                    fprintf(fid, ',%.6f', metrics);
                    fprintf(fid, '\n');
                end

                % Save 9 sample images
                figure('visible','off');
                colormap(fibroclr);
                for k = 1:9
                    subplot(3,3,k);
                    imagesc(example_patterns{k});
                    axis equal off;
                    title(sprintf('%s #%d', type_name, k));
                end
                img_name = sprintf('examples_%s_%s.png', type_name, mode);
                print(img_name, '-dpng');
                close;
            end
        end
    end
end

fclose(fid);
disp('Simulation complete. All metrics saved.');
