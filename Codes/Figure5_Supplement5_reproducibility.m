%% FIGURE 5 + SUPPLEMENT 5 REPRODUCIBILITY ANALYSIS
% =========================================================================
% Reduced, self-contained MATLAB analysis containing only the computations
% used for the data-derived panels in Figure 5 and Supplement 5.
%
% FIGURE 5
%   B-D  Example XY trajectories for WT, AFD-specific eat-4 cKO, and
%        AFD-specific unc-17 cKO.
%   B'   WT plate-level mean position during minutes 10-15, with WT-defined
%        mean +/-1, +/-2, +/-3 SD landmarks.
%   E    Mean temperature position over the 15-min assay (plate hierarchy).
%   F-H  Occupancy heatmaps relative to the six WT-defined landmarks.
%
% SUPPLEMENT 5
%   A-C  Temperature-position histograms for 0-3, 6-9, and 12-15 min.
%   D-I  Plate-level occupancy relative to each WT-defined landmark over time.
%
% Figure 5A/A' are schematics and are not generated from tracking data.
%
% IMPORTANT STATISTICAL UNIT
% --------------------------
% One assay plate is the independent biological replicate. Individual
% trajectories are nested within plates and are not treated as independent
% biological replicates for inferential statistics.
%
% NOTE ABOUT FIGURE 5E
% --------------------
% The supplied long analysis script generates the plate-based temperature
% trajectories but does not contain the statistical calculation that produced
% the single p-value printed on the assembled Figure 5E. To keep this reduced
% script reproducible rather than hard-coding that number, the script computes
% a three-group Kruskal-Wallis test on the plate-level mean temperature during
% minutes 10-15 and annotates that value on the trajectory plot. Verify that
% this is the intended manuscript test before repository release.
%
% MATLAB dependencies for statistics:
%   Statistics and Machine Learning Toolbox (ranksum, kruskalwallis)
% =========================================================================

clc
clear
close all

%% ========================================================================
% USER-EDITABLE PATHS
% ========================================================================

data_dir = ...
    "/Users/andreacuentas/Library/CloudStorage/OneDrive-YaleUniversity/" + ...
    "matlab_behavior_analysis/Steep_15_25_15min_Xcorrect_Tc15_new";

output_dir = fullfile(data_dir,"Figure5_Supplement5_reproducibility_output");

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% ========================================================================
% SETTINGS
% ========================================================================

% Recording / assay.
fps = 2;
assay_duration_min = 15;
nFramesRequired = assay_duration_min*60*fps;

% Thermal calibration.
x_min_px = 900;
x_max_px = 1800;
T_min_C = 15;
T_max_C = 25;

px_to_C = @(x) ...
    T_min_C + (x-x_min_px).*(T_max_C-T_min_C)./(x_max_px-x_min_px);

% Genotype order and plot colors used in the manuscript.
genotype_names = ["wild type","eat-4 cKO","unc-17 cKO"];
plot_colors = [ ...
    0.00 0.00 0.00; ...  % WT
    0.25 0.65 0.38; ...  % eat-4
    0.78 0.25 0.68];     % unc-17

% WT landmark reference window.
reference_start_min = 10;
reference_end_min = 15;
referenceFrames = ...
    round(reference_start_min*60*fps)+1 : ...
    round(reference_end_min*60*fps);

% A trajectory must have valid X measurements for >=50% of the specified
% interval to contribute to the WT landmark/reference-position calculation.
reference_min_valid_fraction = 0.50;

% Full-assay occupancy analysis.
behaviorFrames = 1:nFramesRequired;
behavior_time_min = ((0:nFramesRequired-1)./fps)./60;

% Preserve the validity rule used in the supplied landmark analysis:
% >=50% valid tracking across the assay AND a valid initial X value during
% the first 10 s.
behavior_min_valid_fraction = 0.50;
initial_valid_window_s = 10;
initialFrames = 1:min(nFramesRequired,round(initial_valid_window_s*fps));

% Heatmaps are sampled at 1-min resolution from plate-level occupancy curves.
heatmap_minutes = 0:1:15;
heatmap_frames = round(heatmap_minutes*60*fps)+1;
heatmap_frames = min(heatmap_frames,nFramesRequired);

% Supplement 5 temperature histograms.
% The final supplement displays only these three 3-min windows.
hist_windows_min = [0 3; 6 9; 12 15];
hist_window_labels = ["0-3 min","6-9 min","12-15 min"];
hist_edges_C = T_min_C:0.25:T_max_C;

% Figure 5B-D trajectory examples.
n_example_trajectories = 10;
trajectory_line_width = 1.0;

%% ========================================================================
% LOAD THE THREE DATASETS
% ========================================================================

[X_N2,Y_N2,plate_N2] = local_load_xy_and_plate( ...
    fullfile(data_dir,"X_N2.mat"), ...
    fullfile(data_dir,"Y_N2.mat"), ...
    fullfile(data_dir,"plate_N2.mat"), ...
    "X_N2","Y_N2", ...
    ["plate_N2","plate_id_N2","PLATE_ID_N2","PLATE_ID", ...
     "plate_id","Plate_ID","PlateID","plateID"]);

[X_eat4,Y_eat4,plate_eat4] = local_load_xy_and_plate( ...
    fullfile(data_dir,"X_eat4.mat"), ...
    fullfile(data_dir,"Y_eat4.mat"), ...
    fullfile(data_dir,"plate_eat4.mat"), ...
    "X_eat4","Y_eat4", ...
    ["plate_eat4","plate_id_eat4","PLATE_ID_eat4","PLATE_ID", ...
     "plate_id","Plate_ID","PlateID","plateID"]);

[X_unc17,Y_unc17,plate_unc17] = local_load_xy_and_plate( ...
    fullfile(data_dir,"X_unc17.mat"), ...
    fullfile(data_dir,"Y_unc17.mat"), ...
    fullfile(data_dir,"plate_unc17.mat"), ...
    "X_unc17","Y_unc17", ...
    ["plate_unc17","plate_id_unc17","PLATE_ID_unc17","PLATE_ID", ...
     "plate_id","Plate_ID","PlateID","plateID"]);

X_by_genotype = {X_N2,X_eat4,X_unc17};
Y_by_genotype = {Y_N2,Y_eat4,Y_unc17};
plate_by_genotype = {string(plate_N2(:))',string(plate_eat4(:))', ...
    string(plate_unc17(:))'};

for g = 1:3
    Xg = X_by_genotype{g};
    Yg = Y_by_genotype{g};

    if size(Xg,1) < nFramesRequired || size(Yg,1) < nFramesRequired
        error('%s contains fewer than %d frames required for 15 min.', ...
            genotype_names(g),nFramesRequired);
    end

    X_by_genotype{g} = Xg(1:nFramesRequired,:);
    Y_by_genotype{g} = Yg(1:nFramesRequired,:);

    if ~isequal(size(X_by_genotype{g}),size(Y_by_genotype{g}))
        error('%s: X and Y matrices do not have matching dimensions.', ...
            genotype_names(g));
    end

    if numel(plate_by_genotype{g}) ~= size(X_by_genotype{g},2)
        error('%s: plate vector must contain one entry per trajectory column.', ...
            genotype_names(g));
    end
end

%% ========================================================================
% DEFINE THE WT REFERENCE DISTRIBUTION AND SD LANDMARKS
% ========================================================================
% Hierarchy:
%   frame -> animal mean X during min 10-15
%         -> mean of animal means within each independent WT plate
%         -> WT mean and SD across plate means

[WT_plate_ID,WT_plate_N,WT_plate_mean_X_px,WT_plate_mean_C] = ...
    local_compute_plate_mean_position_window( ...
    X_by_genotype{1},plate_by_genotype{1},referenceFrames, ...
    reference_min_valid_fraction,x_min_px,x_max_px,T_min_C,T_max_C);

WT_mean_px = mean(WT_plate_mean_X_px,'omitnan');
WT_SD_px = std(WT_plate_mean_X_px,0,'omitnan');
WT_mean_C = px_to_C(WT_mean_px);

landmark_px = [ ...
    WT_mean_px-3*WT_SD_px, ...
    WT_mean_px-2*WT_SD_px, ...
    WT_mean_px-1*WT_SD_px, ...
    WT_mean_px+1*WT_SD_px, ...
    WT_mean_px+2*WT_SD_px, ...
    WT_mean_px+3*WT_SD_px];

landmark_C = px_to_C(landmark_px);

landmark_labels = ["-3 SD","-2 SD","-1 SD","+1 SD","+2 SD","+3 SD"];

% Occupancy side:
%   negative SD landmarks -> fraction currently colder than/before landmark
%   positive SD landmarks -> fraction currently warmer than/after landmark
landmark_direction = [-1 -1 -1 +1 +1 +1];

WT_landmarks = table( ...
    landmark_labels(:),landmark_px(:),landmark_C(:),landmark_direction(:), ...
    'VariableNames',{'Landmark','X_px','Temperature_C','Direction'});

writetable(WT_landmarks,fullfile(output_dir,"WT_defined_landmarks.csv"));

WT_plate_summary = table( ...
    WT_plate_ID(:),WT_plate_N(:),WT_plate_mean_X_px(:),WT_plate_mean_C(:), ...
    'VariableNames',{'Plate','N_animals','Mean_X_min10_15_px', ...
    'Mean_T_min10_15_C'});

writetable(WT_plate_summary, ...
    fullfile(output_dir,"Figure5Bprime_WT_plate_positions_10_15min.csv"));

%% ========================================================================
% COMPUTE THE SAME 10-15 MIN PLATE POSITION ENDPOINT FOR ALL GENOTYPES
% ========================================================================

plate_position_10_15 = table();

for g = 1:3
    [pID,pN,pX,pC] = local_compute_plate_mean_position_window( ...
        X_by_genotype{g},plate_by_genotype{g},referenceFrames, ...
        reference_min_valid_fraction,x_min_px,x_max_px,T_min_C,T_max_C);

    T = table( ...
        repmat(genotype_names(g),numel(pID),1),pID(:),pN(:),pX(:),pC(:), ...
        'VariableNames',{'Genotype','Plate','N_animals', ...
        'Mean_X_min10_15_px','Mean_T_min10_15_C'});

    plate_position_10_15 = local_append_table(plate_position_10_15,T);
end

writetable(plate_position_10_15, ...
    fullfile(output_dir,"plate_position_10_15_all_genotypes.csv"));

%% ========================================================================
% FIGURE 5B': WT PLATE-LEVEL POSITION DURING MINUTES 10-15
% ========================================================================

fig = figure('Color','w','Position',[120 120 430 430], ...
    'Name','Figure5Bprime_WT_plate_position_10_15');
ax = axes(fig);
hold(ax,'on');

vals = WT_plate_mean_C(:);
rng(31);
jitter = (rand(size(vals))-0.5)*0.16;

scatter(ax,1+jitter,vals,45,[0.65 0.65 0.65],'filled', ...
    'MarkerFaceAlpha',0.75);

plot(ax,[0.78 1.22],[mean(vals,'omitnan') mean(vals,'omitnan')], ...
    '-','Color',[0 0 0],'LineWidth',2.5);

for k = 1:6
    yline(ax,landmark_C(k),'--',sprintf('%s',landmark_labels(k)), ...
        'Color',[0.55 0.55 0.55],'LineWidth',0.9, ...
        'LabelHorizontalAlignment','right');
end

xlim(ax,[0.55 1.45]);
ylim(ax,[15 19]);
xticks(ax,1);
xticklabels(ax,{'Wild type'});
ylabel(ax,'Mean position during minutes 10-15 (°C)');
title(ax,'Wild type','FontWeight','normal');
box(ax,'off');
set(ax,'TickDir','out');

local_save_pdf(fig,output_dir,"Figure5Bprime_WT_plate_position_10_15.pdf");

%% ========================================================================
% FIGURE 5B-D: REPRESENTATIVE ABSOLUTE XY TRAJECTORIES
% ========================================================================
% The supplied trajectory script uses the first valid trajectories in column
% order when selection_mode="all". We preserve that deterministic rule and
% save the exact selected track columns. Within each genotype, the displayed
% tracks are color-coded with a shade gradient so individual trajectories are
% easier to distinguish.

selected_trajectory_table = table();

% Determine a shared Y range across the selected trajectories after selection.
selected_by_genotype = cell(3,1);
all_selected_Y = [];

for g = 1:3
    joint_valid = isfinite(X_by_genotype{g}) & isfinite(Y_by_genotype{g});
    valid_tracks = find(sum(joint_valid,1)>=2);
    nKeep = min(n_example_trajectories,numel(valid_tracks));
    selected_by_genotype{g} = valid_tracks(1:nKeep);

    YY = Y_by_genotype{g}(:,selected_by_genotype{g});
    all_selected_Y = [all_selected_Y; YY(isfinite(YY))]; %#ok<AGROW>

    selected_trajectory_table = local_append_table( ...
        selected_trajectory_table, ...
        table(repmat(genotype_names(g),nKeep,1), ...
        selected_by_genotype{g}(:), ...
        'VariableNames',{'Genotype','TrackColumn'}));
end

if isempty(all_selected_Y)
    shared_Y_limits = [0 1];
else
    ylo = min(all_selected_Y);
    yhi = max(all_selected_Y);
    pad = max(20,0.03*(yhi-ylo));
    shared_Y_limits = [ylo-pad yhi+pad];
end

trajectory_filenames = [ ...
    "Figure5B_WT_example_trajectories.pdf", ...
    "Figure5C_eat4_example_trajectories.pdf", ...
    "Figure5D_unc17_example_trajectories.pdf"];

for g = 1:3
    fig = local_plot_absolute_trajectories( ...
        X_by_genotype{g},Y_by_genotype{g},selected_by_genotype{g}, ...
        px_to_C,plot_colors(g,:),genotype_names(g), ...
        landmark_C(3),landmark_C(4), ...
        shared_Y_limits,trajectory_line_width);

    local_save_pdf(fig,output_dir,trajectory_filenames(g));
end

writetable(selected_trajectory_table, ...
    fullfile(output_dir,"Figure5BCD_selected_trajectory_columns.csv"));

%% ========================================================================
% PLATE-LEVEL LANDMARK OCCUPANCY ACROSS THE FULL ASSAY
% ========================================================================
% For each plate and time point:
%
% negative landmarks:
%   occupancy = N animals with X <= landmark / N valid animals at that time
%
% positive landmarks:
%   occupancy = N animals with X >= landmark / N valid animals at that time
%
% The plate-level occupancy curve is the fundamental measurement. Genotype
% curves/heatmaps are averages across independent plates.

% Do not preallocate with an empty-field struct. MATLAB requires every
% element of a struct array to have exactly the same fields. The occupancy
% function returns a populated struct, so initialize the array from the
% first returned result and then overwrite subsequent elements.
occupancy_results = [];

for g = 1:3
    this_result = local_compute_plate_occupancy( ...
        X_by_genotype{g},plate_by_genotype{g}, ...
        behaviorFrames,landmark_px,landmark_direction, ...
        behavior_min_valid_fraction,initialFrames,behavior_time_min);

    this_result.Genotype = genotype_names(g);

    if g == 1
        occupancy_results = repmat(this_result,3,1);
    else
        occupancy_results(g) = this_result;
    end
end

%% ========================================================================
% BUILD PLATE-LEVEL OCCUPANCY-AUC SOURCE TABLE
% ========================================================================

occupancy_AUC_table = table();

for g = 1:3
    R = occupancy_results(g);

    for p = 1:numel(R.plateIDs)
        for k = 1:6
            row = table( ...
                genotype_names(g),R.plateIDs(p),landmark_labels(k), ...
                landmark_C(k),R.plateValidN(p),R.plateOccupancyAUC(p,k), ...
                'VariableNames',{'Genotype','Plate','Landmark', ...
                'Landmark_C','N_valid_tracks','Occupancy_AUC_fraction_min'});

            occupancy_AUC_table = local_append_table(occupancy_AUC_table,row);
        end
    end
end

writetable(occupancy_AUC_table, ...
    fullfile(output_dir,"Supplement5_occupancy_AUC_by_plate.csv"));

%% ========================================================================
% SUPPLEMENT 5D-I STATISTICS
% ========================================================================
% Planned comparisons at each landmark:
%   WT vs eat-4 cKO
%   WT vs unc-17 cKO
%
% Wilcoxon rank-sum at the plate level, with Holm correction across those
% two WT-referenced comparisons separately within each landmark.

occupancy_AUC_stats = local_plate_metric_stats( ...
    occupancy_AUC_table,'Occupancy_AUC_fraction_min', ...
    genotype_names,landmark_labels);

writetable(occupancy_AUC_stats, ...
    fullfile(output_dir,"Supplement5_occupancy_AUC_statistics.csv"));

%% ========================================================================
% FIGURE 5F-H: OCCUPANCY HEATMAPS
% ========================================================================

heatmap_data = nan(6,numel(heatmap_minutes),3);

for g = 1:3
    for k = 1:6
        M = occupancy_results(g).plateOccupancy(:,:,k);
        mu = mean(M,1,'omitnan');
        heatmap_data(k,:,g) = mu(heatmap_frames);
    end
end

% Export exact heatmap source data in long format.
HeatGenotype = strings(0,1);
HeatLandmark = strings(0,1);
HeatLandmark_C = zeros(0,1);
HeatTime_min = zeros(0,1);
HeatOccupancy = zeros(0,1);

for g = 1:3
    for k = 1:6
        for t = 1:numel(heatmap_minutes)
            HeatGenotype(end+1,1) = genotype_names(g); %#ok<SAGROW>
            HeatLandmark(end+1,1) = landmark_labels(k); %#ok<SAGROW>
            HeatLandmark_C(end+1,1) = landmark_C(k); %#ok<SAGROW>
            HeatTime_min(end+1,1) = heatmap_minutes(t); %#ok<SAGROW>
            HeatOccupancy(end+1,1) = heatmap_data(k,t,g); %#ok<SAGROW>
        end
    end
end

heatmap_source = table( ...
    HeatGenotype,HeatLandmark,HeatLandmark_C,HeatTime_min,HeatOccupancy, ...
    'VariableNames',{'Genotype','Landmark','Landmark_C','Time_min', ...
    'Mean_plate_occupancy'});

writetable(heatmap_source, ...
    fullfile(output_dir,"Figure5FGH_occupancy_heatmap_source_data.csv"));

fig = figure('Color','w','Position',[80 130 1120 410], ...
    'Name','Figure5FGH_occupancy_heatmaps');
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

for g = 1:3
    ax = nexttile(tl,g);

    % Insert a blank row between -1 SD and +1 SD to display the WT
    % stabilization region separating colder- and warmer-side landmarks.
    Z = [ ...
        heatmap_data(1,:,g); ... % -3
        heatmap_data(2,:,g); ... % -2
        heatmap_data(3,:,g); ... % -1
        nan(1,numel(heatmap_minutes)); ...
        heatmap_data(4,:,g); ... % +1
        heatmap_data(5,:,g); ... % +2
        heatmap_data(6,:,g)];    % +3

    imagesc(ax,heatmap_minutes,1:7,Z,[0 1]);
    set(ax,'YDir','normal');

    yticks(ax,1:7);
    yticklabels(ax,{ ...
        sprintf('-3 SD | %.2f C',landmark_C(1)), ...
        sprintf('-2 SD | %.2f C',landmark_C(2)), ...
        sprintf('-1 SD | %.2f C',landmark_C(3)), ...
        'Stabilization area', ...
        sprintf('+1 SD | %.2f C',landmark_C(4)), ...
        sprintf('+2 SD | %.2f C',landmark_C(5)), ...
        sprintf('+3 SD | %.2f C',landmark_C(6))});

    xlabel(ax,'Time (min)');
    title(ax,genotype_names(g),'FontWeight','bold');

    colormap(ax,parula(256));
    if g==3
        cb=colorbar(ax);
        cb.Label.String='Proportion of animals';
    end

    box(ax,'off');
end

local_save_pdf(fig,output_dir,"Figure5FGH_occupancy_heatmaps.pdf");

%% ========================================================================
% SUPPLEMENT 5A-C: TEMPERATURE-POSITION HISTOGRAMS
% ========================================================================
% These panels are descriptive: valid X samples are pooled within each
% genotype and time window and probability-normalized. They are not used as
% independent observations for inferential statistics.

hist_source = table();
hist_probabilities = cell(size(hist_windows_min,1),3);
hist_centers_C = hist_edges_C(1:end-1)+diff(hist_edges_C)/2;

for w = 1:size(hist_windows_min,1)
    startFrame = round(hist_windows_min(w,1)*60*fps)+1;
    endFrame = round(hist_windows_min(w,2)*60*fps);
    frameIdx = startFrame:endFrame;

    for g = 1:3
        xVals = X_by_genotype{g}(frameIdx,:);
        xVals = xVals(isfinite(xVals));
        tempVals = px_to_C(xVals);

        counts = histcounts(tempVals,hist_edges_C);
        if sum(counts)>0
            probs = counts./sum(counts);
        else
            probs = nan(size(counts));
        end
        hist_probabilities{w,g}=probs;

        tmp = table( ...
            repmat(genotype_names(g),numel(hist_centers_C),1), ...
            repmat(hist_window_labels(w),numel(hist_centers_C),1), ...
            hist_centers_C(:),probs(:), ...
            'VariableNames',{'Genotype','TimeWindow','TemperatureBinCenter_C', ...
            'Probability'});

        hist_source = local_append_table(hist_source,tmp);
    end
end

writetable(hist_source, ...
    fullfile(output_dir,"Supplement5ABC_temperature_histogram_source_data.csv"));

fig = figure('Color','w','Position',[70 70 1080 760], ...
    'Name','Supplement5ABC_temperature_histograms');
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
axs = gobjects(3,3);

for w = 1:3
    rowMax = 0;

    % Determine common vertical range for this row.
    for g = 1:3
        p = hist_probabilities{w,g};
        rowMax = max(rowMax,max(p,[],'omitnan'));
    end
    if ~isfinite(rowMax) || rowMax<=0,rowMax=1;end
    rowMax = rowMax*1.12;

    for g = 1:3
        ax = nexttile(tl,(w-1)*3+g);
        axs(w,g)=ax;

        probs = hist_probabilities{w,g};
        bar(ax,hist_centers_C,probs,1, ...
            'FaceColor',plot_colors(g,:), ...
            'EdgeColor',plot_colors(g,:));

        hold(ax,'on');

        % Outer WT stabilization landmarks shown in the supplement.
        xline(ax,landmark_C(1),'--','-3 SD', ...
            'Color',[0.30 0.35 0.75],'LineWidth',0.9);
        xline(ax,landmark_C(6),'--','+3 SD', ...
            'Color',[0.85 0.35 0.30],'LineWidth',0.9);

        xlim(ax,[15 25]);
        ylim(ax,[0 rowMax]);

        if w==1
            title(ax,genotype_names(g),'FontWeight','bold');
        end

        if g==1
            ylabel(ax,hist_window_labels(w));
        end

        if w==3
            xlabel(ax,'Temperature (°C)');
        else
            ax.XTickLabel=[];
        end

        box(ax,'off');
        set(ax,'TickDir','out');
    end
end

local_save_pdf(fig,output_dir,"Supplement5ABC_temperature_histograms.pdf");

%% ========================================================================
% FIGURE 5E: MEAN TEMPERATURE POSITION ACROSS TIME
% ========================================================================
% Each plate first contributes one mean-temperature trajectory. Genotype
% mean and SEM are then computed across independent plate trajectories.

plate_mean_temp_curves = cell(3,1);
mean_temp_by_genotype = nan(3,nFramesRequired);
sd_temp_by_genotype = nan(3,nFramesRequired);
sem_temp_by_genotype = nan(3,nFramesRequired);
n_temp_plates = nan(3,nFramesRequired);

plate_curve_source = table();

for g = 1:3
    R = occupancy_results(g);
    Xg = X_by_genotype{g};
    plateVec = plate_by_genotype{g};

    plateTemp = nan(numel(R.plateIDs),nFramesRequired);

    for p = 1:numel(R.plateIDs)
        idx = R.validTrack & plateVec==R.plateIDs(p);
        meanXplate = mean(Xg(:,idx),2,'omitnan');
        plateTemp(p,:) = px_to_C(meanXplate)';
    end

    plate_mean_temp_curves{g}=plateTemp;
    mean_temp_by_genotype(g,:)=mean(plateTemp,1,'omitnan');
    sd_temp_by_genotype(g,:)=std(plateTemp,0,1,'omitnan');
    n_temp_plates(g,:)=sum(isfinite(plateTemp),1);
    sem_temp_by_genotype(g,:)=sd_temp_by_genotype(g,:)./sqrt(n_temp_plates(g,:));
    sem_temp_by_genotype(g,n_temp_plates(g,:)<2)=NaN;

    % Long plate-level source table.
    for p = 1:numel(R.plateIDs)
        good = isfinite(plateTemp(p,:));

        tmp = table( ...
            repmat(genotype_names(g),sum(good),1), ...
            repmat(R.plateIDs(p),sum(good),1), ...
            behavior_time_min(good)',plateTemp(p,good)', ...
            'VariableNames',{'Genotype','Plate','Time_min','MeanTemperature_C'});

        plate_curve_source = local_append_table(plate_curve_source,tmp);
    end
end

writetable(plate_curve_source, ...
    fullfile(output_dir,"Figure5E_plate_mean_temperature_trajectories.csv"));

trajectory_summary = table( ...
    behavior_time_min(:), ...
    mean_temp_by_genotype(1,:)',sd_temp_by_genotype(1,:)',sem_temp_by_genotype(1,:)', ...
    mean_temp_by_genotype(2,:)',sd_temp_by_genotype(2,:)',sem_temp_by_genotype(2,:)', ...
    mean_temp_by_genotype(3,:)',sd_temp_by_genotype(3,:)',sem_temp_by_genotype(3,:)', ...
    'VariableNames',{'Time_min', ...
    'WT_mean_C','WT_SD_C','WT_SEM_C', ...
    'eat4_mean_C','eat4_SD_C','eat4_SEM_C', ...
    'unc17_mean_C','unc17_SD_C','unc17_SEM_C'});

writetable(trajectory_summary, ...
    fullfile(output_dir,"Figure5E_temperature_trajectory_summary.csv"));

% Three-group plate-level endpoint test used for the annotation in this
% reduced reproducibility script.
groups = categorical(plate_position_10_15.Genotype,genotype_names,genotype_names);
[p_position_3group,~,kw_stats] = kruskalwallis( ...
    plate_position_10_15.Mean_T_min10_15_C,groups,'off');

position_stats = table( ...
    "Kruskal-Wallis: plate mean temperature, minutes 10-15", ...
    p_position_3group, ...
    'VariableNames',{'Test','P'});

writetable(position_stats, ...
    fullfile(output_dir,"Figure5E_position_statistics.csv"));

fig = figure('Color','w','Position',[100 120 620 460], ...
    'Name','Figure5E_mean_temperature_position_SEM');
ax = axes(fig);
hold(ax,'on');

for g = 1:3
    mu = mean_temp_by_genotype(g,:);
    semv = sem_temp_by_genotype(g,:);

    fill(ax, ...
        [behavior_time_min fliplr(behavior_time_min)], ...
        [mu-semv fliplr(mu+semv)], ...
        plot_colors(g,:), ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');

    plot(ax,behavior_time_min,mu, ...
        'Color',plot_colors(g,:),'LineWidth',2.0, ...
        'DisplayName',genotype_names(g));
end

xlim(ax,[0 15]);
ylim(ax,[15.5 20.2]);
xticks(ax,0:3:15);
xlabel(ax,'Time (minutes)');
ylabel(ax,'Mean position (°C)');
legend(ax,'Location','southwest','Box','off');

text(ax,0.03,0.96,'Shading = \pm SEM across independent plates', ...
    'Units','normalized','VerticalAlignment','top','FontSize',9);

text(ax,0.98,0.08, ...
    sprintf('Kruskal-Wallis (10-15 min)\\np = %.4g',p_position_3group), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom');

box(ax,'off');
set(ax,'TickDir','out');

local_save_pdf(fig,output_dir,"Figure5E_mean_temperature_position_SEM.pdf");

%% ========================================================================
% SUPPLEMENT 5D-I: OCCUPANCY CURVES AT EACH WT LANDMARK
% ========================================================================
% Panel order in the assembled supplement:
%   top row: +3 SD, +2 SD, +1 SD
%   bottom:  -1 SD, -2 SD, -3 SD

panel_order = [6 5 4 3 2 1];

fig = figure('Color','w','Position',[70 80 1040 650], ...
    'Name','Supplement5DI_landmark_occupancy_curves');
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

for ii = 1:6
    k = panel_order(ii);
    ax = nexttile(tl,ii);
    hold(ax,'on');

    for g = 1:3
        M = occupancy_results(g).plateOccupancy(:,:,k);
        local_plot_mean_sem(ax,behavior_time_min,M,plot_colors(g,:), ...
            genotype_names(g));
    end

    xlim(ax,[0 15]);
    ylim(ax,[0 1]);
    xticks(ax,0:5:15);

    if k>=4
        ylabel(ax,'Fraction currently after landmark');
    else
        ylabel(ax,'Fraction currently before landmark');
    end

    xlabel(ax,'Time (min)');
    title(ax,sprintf('%s | %.2f C',landmark_labels(k),landmark_C(k)), ...
        'FontWeight','bold');

    % Annotate significant Holm-adjusted WT comparisons. The statistic is
    % the plate-level AUC of the complete occupancy curve.
    rows = occupancy_AUC_stats.Landmark==landmark_labels(k);
    S = occupancy_AUC_stats(rows,:);

    yText = 0.10;
    for rr = 1:height(S)
        if isfinite(S.P_Holm(rr)) && S.P_Holm(rr)<0.05
            comp = string(S.Comparison(rr));

            if contains(comp,"eat-4")
                txtColor = plot_colors(2,:);
            elseif contains(comp,"unc-17")
                txtColor = plot_colors(3,:);
            else
                txtColor = [0 0 0];
            end

            text(ax,0.98,yText, ...
                sprintf('%s  p_{Holm}=%.3g',comp,S.P_Holm(rr)), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'Color',txtColor,'FontSize',8,'Interpreter','tex');

            yText = yText+0.07;
        end
    end

    if ii==1
        legend(ax,'Location','best','Box','off');
    end

    box(ax,'off');
    set(ax,'TickDir','out');
end

local_save_pdf(fig,output_dir,"Supplement5DI_landmark_occupancy_curves.pdf");

%% ========================================================================
% SAVE A COMPACT MAT FILE WITH THE MANUSCRIPT-FACING ANALYSIS OBJECTS
% ========================================================================

save(fullfile(output_dir,"Figure5_Supplement5_reproducibility_results.mat"), ...
    'WT_landmarks','WT_plate_summary','plate_position_10_15', ...
    'selected_trajectory_table','occupancy_results','occupancy_AUC_table', ...
    'occupancy_AUC_stats','heatmap_data','heatmap_source','hist_source', ...
    'plate_mean_temp_curves','trajectory_summary','position_stats', ...
    '-v7.3');

fprintf('\n============================================================\n');
fprintf('FIGURE 5 / SUPPLEMENT 5 REPRODUCIBILITY ANALYSIS COMPLETE\n');
fprintf('Outputs saved to:\n  %s\n',output_dir);
fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [X,Y,plateID] = local_load_xy_and_plate( ...
    xFile,yFile,plateFile,xVar,yVar,plateCandidates)
% Load X/Y matrices and one plate-ID value per trajectory column.

    SX = load(xFile);
    SY = load(yFile);
    SP = load(plateFile);

    if isfield(SX,xVar)
        X = double(SX.(xVar));
    else
        X = local_first_numeric_matrix(SX);
    end

    if isfield(SY,yVar)
        Y = double(SY.(yVar));
    else
        Y = local_first_numeric_matrix(SY);
    end

    if ~isequal(size(X),size(Y))
        if isequal(size(X),fliplr(size(Y)))
            Y=Y';
        else
            error('X and Y matrices do not have matching dimensions.');
        end
    end

    plateID = [];

    for k=1:numel(plateCandidates)
        nm = char(plateCandidates(k));
        if isfield(SP,nm)
            v = SP.(nm);
            if numel(v)==size(X,2)
                plateID = v;
                break
            end
        end
    end

    if isempty(plateID)
        f = fieldnames(SP);
        for k=1:numel(f)
            v=SP.(f{k});
            if (isnumeric(v) || isstring(v) || iscellstr(v) || iscategorical(v)) ...
                    && numel(v)==size(X,2)
                plateID=v;
                break
            end
        end
    end

    if isempty(plateID)
        error('Could not find a plate-ID vector with one entry per trajectory.');
    end

    plateID = string(plateID(:));
end


function M = local_first_numeric_matrix(S)
    f=fieldnames(S);
    for k=1:numel(f)
        v=S.(f{k});
        if isnumeric(v) && ismatrix(v) && ~isscalar(v)
            M=double(v);
            return
        end
    end
    error('No suitable numeric matrix was found in MAT file.');
end


function [plateIDs,plateNAnimals,plateMeanXpx,plateMeanC] = ...
    local_compute_plate_mean_position_window( ...
    X,plateID,analysisFrames,minValidFraction, ...
    xMinPx,xMaxPx,TminC,TmaxC)
% One value per independent plate for a position-analysis window.
%
% frame -> animal mean -> plate mean

    if size(X,1)<max(analysisFrames)
        error('Position analysis requires more frames than are available.');
    end

    plateID=string(plateID(:))';

    animalMeanX=mean(X(analysisFrames,:),1,'omitnan');
    fractionValid=mean(isfinite(X(analysisFrames,:)),1);

    validAnimal=fractionValid>=minValidFraction & isfinite(animalMeanX);

    plateIDs=unique(plateID(validAnimal),'stable')';
    plateNAnimals=nan(numel(plateIDs),1);
    plateMeanXpx=nan(numel(plateIDs),1);

    for p=1:numel(plateIDs)
        idx=validAnimal & plateID==plateIDs(p);
        plateNAnimals(p)=sum(idx);
        plateMeanXpx(p)=mean(animalMeanX(idx),'omitnan');
    end

    good=isfinite(plateMeanXpx);
    plateIDs=plateIDs(good);
    plateNAnimals=plateNAnimals(good);
    plateMeanXpx=plateMeanXpx(good);

    plateMeanC = ...
        TminC + (plateMeanXpx-xMinPx).*(TmaxC-TminC)./(xMaxPx-xMinPx);
end


function R = local_compute_plate_occupancy( ...
    X,plateID,assayFrames,landmarkPx,landmarkDirection, ...
    minValidFraction,initialFrames,timeMin)
% Compute landmark occupancy using one independent curve per plate.
%
% Valid trajectory rule preserved from the supplied landmark analysis:
%   >= minValidFraction valid X across the assay
%   AND a finite X value within the first initialFrames.
%
% negative landmark: occupancy = fraction currently X <= landmark
% positive landmark: occupancy = fraction currently X >= landmark

    X=X(assayFrames,:);
    nFrames=size(X,1);
    nAnimals=size(X,2);
    nLandmarks=numel(landmarkPx);

    plateID=string(plateID(:))';
    if numel(plateID)~=nAnimals
        error('Plate-ID length does not match number of trajectory columns.');
    end

    validFraction=mean(isfinite(X),1);

    initialX=nan(1,nAnimals);
    for j=1:nAnimals
        q=find(isfinite(X(initialFrames,j)),1,'first');
        if ~isempty(q)
            initialX(j)=X(initialFrames(q),j);
        end
    end

    validTrack=validFraction>=minValidFraction & isfinite(initialX);
    plateIDs=unique(plateID(validTrack),'stable');

    nPlates=numel(plateIDs);
    plateValidN=zeros(nPlates,1);
    plateOccupancy=nan(nPlates,nFrames,nLandmarks);
    plateOccupancyAUC=nan(nPlates,nLandmarks);

    for p=1:nPlates
        pIdx=validTrack & plateID==plateIDs(p);
        plateValidN(p)=sum(pIdx);
        Xp=X(:,pIdx);

        for k=1:nLandmarks
            if landmarkDirection(k)>0
                beyond=Xp>=landmarkPx(k);
            else
                beyond=Xp<=landmarkPx(k);
            end

            validNow=isfinite(Xp);
            numerator=sum(beyond & validNow,2);
            denominator=sum(validNow,2);

            occ=numerator./denominator;
            occ(denominator==0)=NaN;

            plateOccupancy(p,:,k)=occ';

            good=isfinite(occ) & isfinite(timeMin(:));
            if sum(good)>=2
                plateOccupancyAUC(p,k)=trapz(timeMin(good),occ(good)');
            end
        end
    end

    R=struct();
    R.plateIDs=plateIDs(:);
    R.validTrack=validTrack;
    R.validFraction=validFraction;
    R.initialX=initialX;
    R.plateValidN=plateValidN;
    R.plateOccupancy=plateOccupancy;
    R.plateOccupancyAUC=plateOccupancyAUC;
end


function fig = local_plot_absolute_trajectories( ...
    X,Y,selectedTracks,pxToC,lineColor,genotypeName, ...
    minus1C,plus1C,yLimits,lineWidth)
% Plot raw absolute XY trajectories with X converted to temperature.
% Individual tracks are displayed with genotype-specific shade gradients.

    fig=figure('Color','w','Position',[120 120 520 510], ...
        'Name',"Trajectory_"+genotypeName);
    ax=axes(fig);
    hold(ax,'on');

    nTracks=numel(selectedTracks);
    trackColors=local_make_trajectory_shade_palette(genotypeName,nTracks,lineColor);

    for ii=1:nTracks
        j=selectedTracks(ii);

        good=isfinite(X(:,j)) & isfinite(Y(:,j));
        xx=pxToC(X(good,j));
        yy=Y(good,j);

        if isempty(xx),continue;end

        thisColor=trackColors(ii,:);

        plot(ax,xx,yy,'Color',thisColor,'LineWidth',lineWidth, ...
            'HandleVisibility','off');

        scatter(ax,xx(1),yy(1),56,thisColor,'filled', ...
            'HandleVisibility','off');

        scatter(ax,xx(end),yy(end),64,thisColor,'o','LineWidth',1.8, ...
            'HandleVisibility','off');
    end

    % WT stabilization-area boundaries.
    xline(ax,minus1C,'--','-1 SD','Color',[0.35 0.35 0.35], ...
        'LineWidth',1.0,'LabelVerticalAlignment','top');
    xline(ax,plus1C,'--','+1 SD','Color',[0.35 0.35 0.35], ...
        'LineWidth',1.0,'LabelVerticalAlignment','top');

    xlim(ax,[15 25]);
    ylim(ax,yLimits);
    xlabel(ax,'Temperature (°C)');
    ylabel(ax,'Y position (pixels)');
    title(ax,genotypeName,'FontWeight','bold');

    % Small legend showing marker meaning.
    scatter(ax,nan,nan,56,[0 0 0],'filled','DisplayName','Start position');
    scatter(ax,nan,nan,64,[0 0 0],'o','LineWidth',1.8, ...
        'DisplayName','Final position');
    legend(ax,'Location','northeast','Box','off');

    box(ax,'off');
    set(ax,'TickDir','out');
end


function palette = local_make_trajectory_shade_palette(genotypeName,nTracks,baseColor)
% Build a genotype-specific shade palette so the example trajectories can
% be visually distinguished within each genotype panel.
%
% WT      : black to light gray
% eat-4   : dark green to light green
% unc-17  : light magenta to dark magenta

    if nTracks<=0
        palette=zeros(0,3);
        return
    end

    gname=lower(strtrim(char(genotypeName)));

    if contains(gname,'wild') || contains(gname,'n2')
        c0=[0.00 0.00 0.00];
        c1=[0.72 0.72 0.72];
    elseif contains(gname,'eat')
        c0=[0.05 0.35 0.12];
        c1=[0.68 0.86 0.68];
    elseif contains(gname,'unc')
        c0=[0.90 0.70 0.88];
        c1=[0.55 0.08 0.45];
    else
        c0=max(baseColor*0.55,0);
        c1=1-(1-baseColor)*0.35;
    end

    t=linspace(0,1,nTracks)';
    palette=(1-t).*c0 + t.*c1;
end


function local_plot_mean_sem(ax,x,M,lineColor,displayName)
% Mean +/- SEM across independent plate curves.

    mu=mean(M,1,'omitnan');
    n=sum(isfinite(M),1);
    sem=std(M,0,1,'omitnan')./sqrt(n);
    sem(n<2)=NaN;

    lo=mu-sem;
    hi=mu+sem;

    fill(ax,[x fliplr(x)],[lo fliplr(hi)],lineColor, ...
        'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');

    plot(ax,x,mu,'Color',lineColor,'LineWidth',2, ...
        'DisplayName',displayName);
end


function statsTable = local_plate_metric_stats( ...
    summaryTable,metricName,genotypeNames,landmarkLabels)
% Planned WT-referenced Wilcoxon rank-sum tests.
%
% Within each landmark:
%   WT vs eat-4
%   WT vs unc-17
%
% Holm correction is applied across these two comparisons.

    Landmark = strings(0,1);
    Metric = strings(0,1);
    Comparison = strings(0,1);
    N_WT = zeros(0,1);
    N_mutant = zeros(0,1);
    P_raw = zeros(0,1);
    P_Holm = zeros(0,1);

    for k=1:numel(landmarkLabels)
        wt = summaryTable.(metricName)( ...
            summaryTable.Landmark==landmarkLabels(k) & ...
            summaryTable.Genotype==genotypeNames(1));
        wt=wt(isfinite(wt));

        raw=nan(2,1);
        nMut=zeros(2,1);

        for q=2:3
            mut=summaryTable.(metricName)( ...
                summaryTable.Landmark==landmarkLabels(k) & ...
                summaryTable.Genotype==genotypeNames(q));
            mut=mut(isfinite(mut));
            nMut(q-1)=numel(mut);

            if numel(wt)>=2 && numel(mut)>=2
                raw(q-1)=ranksum(wt,mut);
            end
        end

        adj=local_holm_adjust(raw);

        for q=2:3
            Landmark(end+1,1)=landmarkLabels(k); %#ok<AGROW>
            Metric(end+1,1)=string(metricName); %#ok<AGROW>
            Comparison(end+1,1)=genotypeNames(1)+" vs "+genotypeNames(q); %#ok<AGROW>
            N_WT(end+1,1)=numel(wt); %#ok<AGROW>
            N_mutant(end+1,1)=nMut(q-1); %#ok<AGROW>
            P_raw(end+1,1)=raw(q-1); %#ok<AGROW>
            P_Holm(end+1,1)=adj(q-1); %#ok<AGROW>
        end
    end

    statsTable=table( ...
        Landmark,Metric,Comparison,N_WT,N_mutant,P_raw,P_Holm, ...
        'VariableNames',{'Landmark','Metric','Comparison', ...
        'N_WT_plates','N_mutant_plates','P_raw','P_Holm'});
end


function adjustedP = local_holm_adjust(rawP)
% Holm step-down family-wise error correction.

    rawP=rawP(:);
    adjustedP=nan(size(rawP));

    good=isfinite(rawP);
    p=rawP(good);

    if isempty(p),return;end

    [ps,ord]=sort(p);
    m=numel(ps);

    adjSorted=nan(m,1);
    running=0;

    for i=1:m
        this=(m-i+1)*ps(i);
        running=max(running,this);
        adjSorted(i)=min(running,1);
    end

    tmp=nan(m,1);
    tmp(ord)=adjSorted;

    adjustedP(good)=tmp;
end


function T = local_append_table(T,row)
    if isempty(T) || width(T)==0
        T=row;
    else
        T=[T;row]; %#ok<AGROW>
    end
end


function local_save_pdf(fig,outputDir,fileName)
    exportgraphics(fig,fullfile(outputDir,fileName), ...
        'ContentType','vector','BackgroundColor','white');
end
