function FIGURE_analyseConductivity(pattern_array, regenerate)
% This function takes the input array of microfibrosis patterns and
% evaluates the diffusivity (in the faster principal direction), and the
% anisotropy ratio (in terms of diffusivity)

% Specify the size of the elements (number of pixels)
Nh = 200;

% Lay out colours
plot_colors = 1/255 * [ [112, 122, 140];
                        [012, 142, 207];
                        [209, 188, 050];
                        [224, 126, 160] ];

% List of legend entries
legend_txt = {'Interstitial', 'Compact', 'Diffuse', 'Patchy'};
                    
% List of patterns to plot
plot_patterns = [1,2,3,4];

% Plotting settings for anisotropy value
AR_max = 20;         % Pin larger values to this value
AR_step = 2;         % Separation between ticks on AR axes
AR_base = 3;         % Base AR value before effects of fibrosis

% Specify the base value of longitudinal conductivity
base_longcond = 3;

% Specify filename to use for storage
filename = 'fullpattern_Ddata';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume data will not be regenerated
if nargin < 2
    regenerate = false;
end

% Create data if not found, or asked to regenerate
if ~exist([filename,'.mat'],'file') || regenerate
        
    % Initialise waitbar
    wtbar_obj = waitbar(0,'Working...');
    
    % Loop over each pattern type
    for k = 1:size(pattern_array,1)
                
        % Loop over each individual pattern, converting each into a problem
        % object so that the homogenisation code may be applied
        parfor j = 1:size(pattern_array,2)
            
            % Create a problem array using the current map
            problem = createProblemFromMap(pattern_array{k,j});
            
            % Homogenise this problem (with specified skin map)
            homog_problems{j} = homogeniseFull2DProblem( problem, Nh, Nh, 'linear' );
            
        end
        
        % Read out the dimensions of the homogenised problems
        [Ex,Ey] = size(homog_problems{1}.D_tensor.D_xx);
        
        % Now loop over each homogenised problem and grab out the
        % conductive properties
        for j = 1:size(pattern_array,2)
                        
            % Process each element in the homogenised problem
            c = 0;
            for ii = 2:Ex-1
                for jj = 2:Ey-1
                    
                    % Increment counter
                    c = c + 1;
                    
                    % Grab out the tensor here
                    D = [homog_problems{j}.D_tensor.D_xx(ii,jj), homog_problems{j}.D_tensor.D_xy(ii,jj); homog_problems{j}.D_tensor.D_xy(ii,jj), homog_problems{j}.D_tensor.D_yy(ii,jj)];
                        
                    %%% Store summaries of conductivity
                                       
                    % Grab out principal directions and strengths
                    [V, D_vals] = eig(D);
                    D_vals = diag(D_vals);
                    % Ensure no negative eigenvalues have appeared
                    D_vals(D_vals < 0) = 0;
                    % Shift all eigenvectors to point left-to-right
                    V = V .* sign(V(1,:));
                    % Find angles of each to the horizontal
                    thetas = acos([1,0] * V);
                    % Find the direction closest to horizontal, and how far
                    % it deviates from zero
                    [close_deviation,closest] = min(abs(thetas));
                    % Find the other index (ugly, but this also works in 
                    % the case of a repeated angle)
                    indices = [1,2];
                    indices(closest) = [];
                    distant = indices(1);
   
                    %%% Summary One - Conductivity in fibre direction
                    fibre_conds_here(c) = [1;0]' * D * [1;0];
                    %%% Summary Two - Anisotropy ratio (closest:distant)
                    ARs_here(c) = D_vals(closest) / D_vals(distant);              
                    %%% Summary Three - Angle away from horizontal (of
                    %%% closest principal direction to horizontal)
                    deviations_here(c) = close_deviation;
                    
                end
            end
            
            fibre_conds(:,j) = fibre_conds_here(:);
            ARs(:,j) = ARs_here(:);
            deviations(:,j) = deviations_here(:);
            
        end
        
        % Group together all these values across the different problems
        fibre_cond(k,:) = fibre_conds(:);
        AR(k,:) = ARs(:);
        deviation(k,:) = deviations(:);
        
        waitbar((k-1)/size(pattern_array,1),wtbar_obj,'Working...');
                
    end
    
    % Close the waitbar
    close(wtbar_obj);
    
    % Save data
    save([filename,'.mat'],'fibre_cond','AR','deviation');
    
else
    
    load([filename,'.mat'],'fibre_cond','AR','deviation');
    
end

% Set all anisotropy values above the maximum value to the max value
AR(AR > AR_max) = AR_max;

% Calculate size of data to be plotted
N_ele = size(AR,2);
N_plot = length(plot_patterns) * N_ele;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FIRST PLOT: LONGITUDINAL CONDUCTIVITY VS ANISOTROPY RATIO %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise scatter data (this setup randomises which dots appear on top
% of which dots, preventing one class from partially hiding another)
scatter_data = zeros(N_plot, 2);
scatter_classes = zeros(N_plot, 1);

% Prepare scatter data
c = 0;
for k = plot_patterns
    
    % Increment class counter
    c = c + 1;
    
    % Add the current class to the scatter data
    scatter_data( (c-1)*N_ele+1:c*N_ele, 1 ) = fibre_cond(k,:);
    scatter_data( (c-1)*N_ele+1:c*N_ele, 2 ) = AR(k,:);
    scatter_classes( (c-1)*N_ele+1:c*N_ele ) = c;
    
end

% Randomise ordering of data
I = randperm(N_plot);
scatter_data = scatter_data(I,:);
scatter_classes = scatter_classes(I);

% Initialise figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
hold on;

% Plot the scatterplot
scatter( scatter_data(:,1), scatter_data(:,2), 75, scatter_classes, 'filled' );
colormap(plot_colors(1:length(plot_patterns),:));

% Plot dummy data to create the legend
c = 0;
for k = plot_patterns
    c = c + 1;
    plot(NaN, NaN, '.', 'MarkerSize', 30, 'MarkerEdgeColor', plot_colors(c,:));
end

% Square axis before further adjustments
axis square;

% Add labels
xlabel('Longitudinal Conductivity','FontSize',24);
ylabel('Anisotropy Ratio','FontSize',24);

% Adjust AR axis (starts at no anisotropy (one), ends at given max value)
xlim([0 base_longcond]);
ylim([0 AR_max]);

% Specify ticks manually
c = 0;
% Increment values up to the base
for k = 0:AR_step:AR_base*0.99
    c = c+1;
    ytick_txt{c} = num2str(k);
    ytick_vals(c) = k;
end
% Add in the base value
c = c+1;
ytick_txt{c} = 'Base';
ytick_vals(c) = AR_base;
% Increment values up to the maximum value now
for k = (floor(AR_base/AR_step)+1)*AR_step:AR_step:AR_max
    c = c+1;
    ytick_txt{c} = num2str(k);
    ytick_vals(c) = k;
end
ytick_txt{end} = [ytick_txt{end},'+'];
set(gca,'ytick',ytick_vals);
set(gca,'yticklabel',ytick_txt);

% Thicken axis and increase fontsize
set(gca,'FontSize',24);
set(gca,'LineWidth',2.5);

% Plot guideline for base AR value
plot([-1,base_longcond+1],[AR_base, AR_base],'--','LineWidth',2,'Color',[0.4 0.4 0.4]);

% Add legend
legend({'',legend_txt{plot_patterns}});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SECOND PLOT: LONGITUDINAL CONDUCTIVITY VS PRINCIPAL DIRECTION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise scatter data (this setup randomises which dots appear on top
% of which dots, preventing one class from partially hiding another)
scatter_data = zeros(N_plot, 2);
scatter_classes = zeros(N_plot, 1);

% Prepare scatter data
c = 0;
for k = plot_patterns
    
    % Increment class counter
    c = c + 1;
    
    % Add the current class to the scatter data
    scatter_data( (c-1)*N_ele+1:c*N_ele, 1 ) = fibre_cond(k,:);
    scatter_data( (c-1)*N_ele+1:c*N_ele, 2 ) = deviation(k,:) * 180/pi;
    scatter_classes( (c-1)*N_ele+1:c*N_ele ) = c;
    
end

% Randomise ordering of data
I = randperm(N_plot);
scatter_data = scatter_data(I,:);
scatter_classes = scatter_classes(I);

% Initialise figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
hold on;

% Plot the scatterplot
scatter( scatter_data(:,1), scatter_data(:,2), 75, scatter_classes, 'filled' );
colormap(plot_colors(1:length(plot_patterns),:));

% Plot dummy data to create the legend
c = 0;
for k = plot_patterns
    c = c + 1;
    plot(NaN, NaN, '.', 'MarkerSize', 30, 'MarkerEdgeColor', plot_colors(c,:));
end


% Square axis before further adjustments
axis square;

% Add labels
xlabel('Longitudinal Conductivity','FontSize',24);
ylabel('Deviation from Fibre Direction (^o)','FontSize',24);

% Adjust AR axis (starts at no anisotropy (one), ends at given max value)
xlim([0 base_longcond]);
ylim([0 45]);

% Thicken axis and increase fontsize
set(gca,'FontSize',24);
set(gca,'LineWidth',2.5);

% Add legend
legend({'',legend_txt{plot_patterns}});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% THIRD PLOT: ANISOTROPY RATIO VS PRINCIPAL DIRECTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise scatter data (this setup randomises which dots appear on top
% of which dots, preventing one class from partially hiding another)
scatter_data = zeros(N_plot, 2);
scatter_classes = zeros(N_plot, 1);

% Prepare scatter data
c = 0;
for k = plot_patterns
    
    % Increment class counter
    c = c + 1;
    
    % Add the current class to the scatter data
    scatter_data( (c-1)*N_ele+1:c*N_ele, 1 ) = AR(k,:);
    scatter_data( (c-1)*N_ele+1:c*N_ele, 2 ) = deviation(k,:) * 180/pi;
    scatter_classes( (c-1)*N_ele+1:c*N_ele ) = c;
    
end

% Randomise ordering of data
I = randperm(N_plot);
scatter_data = scatter_data(I,:);
scatter_classes = scatter_classes(I);

% Initialise figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
hold on;

% Plot the scatterplot
scatter( scatter_data(:,1), scatter_data(:,2), 75, scatter_classes, 'filled' );
colormap(plot_colors(1:length(plot_patterns),:));

% Plot dummy data to create the legend
c = 0;
for k = plot_patterns
    c = c + 1;
    plot(NaN, NaN, '.', 'MarkerSize', 30, 'MarkerEdgeColor', plot_colors(c,:));
end


% Square axis before further adjustments
axis square;

% Add labels
xlabel('Anisotropy Ratio','FontSize',24);
ylabel('Deviation from Fibre Direction (^o)','FontSize',24);

% Adjust AR axis (starts at no anisotropy (one), ends at given max value)
xlim([0 AR_max]);
ylim([0 45]);

% Specify ticks manually
c = 0;
for k = 0:AR_step:AR_max
    c = c+1;
    xtick_txt{c} = num2str(k);
end
xtick_txt{end} = [xtick_txt{end},'+'];
set(gca,'xtick',0:AR_step:AR_max);
set(gca,'xticklabel',xtick_txt);

% Thicken axis and increase fontsize
set(gca,'FontSize',24);
set(gca,'LineWidth',2.5);

% Add legend
legend({'',legend_txt{plot_patterns}});