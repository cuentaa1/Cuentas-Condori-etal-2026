%% FIGURE 6 + SUPPLEMENT 6 REPRODUCIBILITY ANALYSIS
% -------------------------------------------------------------------------
% Purpose
% -------
% Minimal analysis script containing only the computations used for the
% data-derived panels in Figure 6 and Supplement 6.
%
% Included analyses:
%
% Figure 6
%   A  Signed thermal migration velocity: WT vs eat-4 cKO
%   B  Signed thermal migration velocity: WT vs unc-17 cKO
%   C  Reorientation rate vs preceding 2-s dT/dt, 19-21 C: WT vs eat-4 cKO
%   D-E  Longest coldward / warmward runs, 19-21 C: WT vs eat-4 cKO
%   F  Reorientation rate vs preceding 2-s dT/dt, 17-19 C: WT vs unc-17 cKO
%   G-H  Longest coldward / warmward runs, 17-19 C: WT vs unc-17 cKO
%
% Supplement 6
%   C-D  True coldward- and warmward-exit rates, 19-21 C: WT vs eat-4 cKO
%   E    Reorientation transition matrices, 19-21 C: WT vs eat-4 cKO
%   F-G  True coldward- and warmward-exit rates, 17-19 C: WT vs unc-17 cKO
%   H    Reorientation transition matrices, 17-19 C: WT vs unc-17 cKO
%
% Supplement 6A and 6B are schematic/cartoon panels and are NOT generated
% from tracking data by this script.
%
% REQUIREMENT
% -----------
% This script uses visulaize_tracks_runs_and_turns_5 from the
% worm-navigation-strategy analysis code used for trajectory preprocessing.
%
% MATLAB release used during development: R2025a.
% -------------------------------------------------------------------------

clc
clear
close all

%% ========================================================================
% USER-EDITABLE PATHS
% ========================================================================
% Point these to your local copies before running.

worm_navigation_toolbox_dir = ...
    "/Users/andreacuentas/Documents/Matlab/worm-navigation-strategy";

data_dir = ...
    "/Users/andreacuentas/Library/CloudStorage/OneDrive-YaleUniversity/" + ...
    "matlab_behavior_analysis/Steep_15_25_15min_Xcorrect_Tc15_new";

output_dir = fullfile(data_dir,"Figure6_Supplement6_reproducibility_output");

addpath(genpath(worm_navigation_toolbox_dir));

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% ========================================================================
% LOAD TRACKING DATA + PLATE IDENTITIES
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

Xdata = {X_N2,X_eat4,X_unc17};
Ydata = {Y_N2,Y_eat4,Y_unc17};
plateIDdata = {string(plate_N2(:)),string(plate_eat4(:)),string(plate_unc17(:))};

strain_names = {'N2','eat-4 cKO','unc-17 cKO'};

% WT/N2 = black throughout.
strain_colors = [ ...
    0.00 0.00 0.00; ...
    0.00 0.32 0.12; ...
    0.48 0.05 0.62];

nStrains = numel(strain_names);

%% ========================================================================
% ANALYSIS SETTINGS USED FOR FIGURE 6 / SUPPLEMENT 6
% ========================================================================

% Assay:
%   X = 900 px -> 15 C
%   X = 1800 px -> 25 C
%   X increases toward warmer temperatures.
pix_range = [900 1800];
temp_range = [15 25];

sampling_rate_raw = 2;        % raw tracking frames/s
sampling_rate_processed = 1;  % processed/smoothed coordinates

minutes_total = 15;
required_frames = minutes_total*60*sampling_rate_raw;

% Whole-gradient thermal-velocity analysis.
temp_edges_C = 15:1:25;
x_edges = pix_range(1) + ...
    (temp_edges_C-temp_range(1)).*diff(pix_range)./diff(temp_range);
x_edges = round(x_edges);

nXbins = numel(x_edges)-1;
temp_bin_labels = cell(nXbins,1);
for b = 1:nXbins
    temp_bin_labels{b} = sprintf('%.0f-%.0f C',temp_edges_C(b),temp_edges_C(b+1));
end

% A temperature bin enters the thermal-velocity LME only when BOTH
% genotypes have at least four independent assay plates represented.
min_plates_per_genotype_temperature_bin = 4;

% Focused, data-guided ROIs used in the manuscript.
eat4_roi_C = [19 21];
unc17_roi_C = [17 19];

% Reorientation episode and sensory-history definitions.
reorientation_cluster_gap_s = 2;
reorientation_cluster_mode = 'first';
preceding_dTdt_window_s = 2;
heading_window_s = 2;
dTdt_edges_C_per_s = -0.03:0.01:0.03;

% A plate x dT/dt bin contributes if it has >=10 s of qualifying exposure.
min_state_exposure_s = 10;

% A dT/dt bin enters the genotype x dT/dt model only if BOTH genotypes
% contribute >=4 independent plates after the exposure filter.
min_plates_per_genotype_dTdt_bin = 4;

% Direction sectors:
%   abs(heading) < 64 deg        -> Coldward
%   64 <= abs(heading) <= 116   -> Transverse
%   abs(heading) > 116 deg       -> Warmward
coldward_halfwidth_deg = 64;
warmward_threshold_deg = 116;

% Complete-run definition used for the trajectory panels.
min_run_seconds = 5;
run_speed_cutoff_mm_s = 0.04;

% Figure 6 D/E/G/H: plot the longest X runs per genotype x direction.
top_directional_runs_to_plot = 30;

% Direction colors used by the standalone trajectory output.
coldward_run_color = [0.24 0.83 0.99];
warmward_run_color = [1.00 0.60 0.27];

%% ========================================================================
% BASIC INPUT CHECKS + PLATE ORDER
% ========================================================================

plateOrder = cell(1,nStrains);

for s = 1:nStrains
    if size(Xdata{s},1) < required_frames || size(Ydata{s},1) < required_frames
        error('%s has fewer than %d raw frames required for 15 min.', ...
            strain_names{s},required_frames);
    end

    if size(Xdata{s},2) ~= size(Ydata{s},2)
        error('%s: X and Y contain different numbers of trajectory columns.', ...
            strain_names{s});
    end

    if numel(plateIDdata{s}) ~= size(Xdata{s},2)
        error('%s: plate-ID vector must contain one value per trajectory column.', ...
            strain_names{s});
    end

    plateOrder{s} = unique(plateIDdata{s},'stable');
end

%% ========================================================================
% PROCESS FULL 0-15 MIN TRAJECTORIES
% ========================================================================
% All manuscript analyses below use the same processed/smoothed trajectories.

coords_by_strain = cell(1,nStrains);

for s = 1:nStrains
    x = Xdata{s}(1:required_frames,:);
    y = Ydata{s}(1:required_frames,:);

    [~,~,~,~,coords] = visulaize_tracks_runs_and_turns_5( ...
        x,y,strain_names{s},strain_colors(s,:),0);

    if numel(coords) ~= numel(plateIDdata{s})
        error(['%s: processed coordinate count (%d) does not match ' ...
               'trajectory/plate-ID count (%d).'], ...
            strain_names{s},numel(coords),numel(plateIDdata{s}));
    end

    coords_by_strain{s} = coords;
end

%% ========================================================================
% FIGURE 6A-B: SIGNED THERMAL MIGRATION VELOCITY
% ========================================================================
% Every consecutive valid processed coordinate pair contributes.
% The pair is assigned to a 1-C temperature bin using midpoint X.
%
% signed thermal velocity = dT/dt
%   negative -> colder / coldward progress
%   positive -> warmer / warmward progress
%
% Interval values are averaged within track x temperature bin.
% Plate-level descriptive points = mean of track-level means within plate/bin.

thermal_track_table = local_extract_signed_thermal_velocity_tracks_only( ...
    coords_by_strain,plateIDdata,strain_names, ...
    x_edges,temp_bin_labels,temp_range,pix_range,sampling_rate_processed);

thermal_track_table = local_sort_temperature_table(thermal_track_table);

thermal_plate_summary = local_track_table_to_plate_curve_summary( ...
    thermal_track_table,'TrackMeanSignedThermalVelocity_C_per_min');

thermal_pairs = [1 2; 1 3];
thermal_tags = ["Figure6A_WT_vs_eat4","Figure6B_WT_vs_unc17"];

thermal_results = struct([]);

for pp = 1:size(thermal_pairs,1)
    compareIdx = thermal_pairs(pp,:);

    [~,interactionT,overallT,binT,modelInfoT,~,~,analysisT] = ...
        local_fit_signed_track_metric_lme( ...
        thermal_track_table,strain_names,compareIdx,temp_bin_labels, ...
        min_plates_per_genotype_temperature_bin,'holm', ...
        'Signed thermal migration velocity');

    pairResult = struct( ...
        'Interaction',interactionT, ...
        'Overall',overallT, ...
        'BinContrasts',binT, ...
        'ModelInfo',modelInfoT);

    fig = local_plot_paper_thermal_velocity_pair( ...
        thermal_plate_summary,pairResult, ...
        strain_names,strain_colors,compareIdx, ...
        min_plates_per_genotype_temperature_bin);

    local_save_clean_figure(fig, ...
        fullfile(output_dir,thermal_tags(pp)+".pdf"));

    prefix = fullfile(output_dir,thermal_tags(pp));
    writetable(binT,prefix+"_bin_statistics.csv");
    writetable(overallT,prefix+"_overall_statistics.csv");
    writetable(modelInfoT,prefix+"_model_info.csv");
    writetable(analysisT,prefix+"_model_input.csv");

    thermal_results(pp).CompareIdx = compareIdx; %#ok<SAGROW>
    thermal_results(pp).Interaction = interactionT;
    thermal_results(pp).Overall = overallT;
    thermal_results(pp).BinContrasts = binT;
end

writetable(thermal_track_table, ...
    fullfile(output_dir,"Figure6AB_signed_thermal_velocity_track_source_data.csv"));
writetable(thermal_plate_summary, ...
    fullfile(output_dir,"Figure6AB_signed_thermal_velocity_plate_source_data.csv"));

%% ========================================================================
% FIGURE 6C-E + SUPPLEMENT 6C-E: EAT-4 ROI, 19-21 C
% ========================================================================

eat4_results = local_run_paper_focused_roi( ...
    coords_by_strain,plateIDdata,plateOrder, ...
    strain_names,strain_colors,[1 2],eat4_roi_C, ...
    sampling_rate_processed,pix_range,temp_range, ...
    preceding_dTdt_window_s,heading_window_s,dTdt_edges_C_per_s, ...
    min_state_exposure_s,min_plates_per_genotype_dTdt_bin, ...
    reorientation_cluster_gap_s,reorientation_cluster_mode, ...
    coldward_halfwidth_deg,warmward_threshold_deg, ...
    min_run_seconds,run_speed_cutoff_mm_s, ...
    top_directional_runs_to_plot,coldward_run_color,warmward_run_color, ...
    output_dir,"eat4_19_21C", ...
    "Figure6C_reorientation_vs_preceding_dTdt_WT_vs_eat4", ...
    "Figure6DE_longest_runs_WT_vs_eat4", ...
    "Supplement6CD_true_exit_rates_WT_vs_eat4", ...
    "Supplement6E_transition_matrices_WT_vs_eat4");

%% ========================================================================
% FIGURE 6F-H + SUPPLEMENT 6F-H: UNC-17 ROI, 17-19 C
% ========================================================================

unc17_results = local_run_paper_focused_roi( ...
    coords_by_strain,plateIDdata,plateOrder, ...
    strain_names,strain_colors,[1 3],unc17_roi_C, ...
    sampling_rate_processed,pix_range,temp_range, ...
    preceding_dTdt_window_s,heading_window_s,dTdt_edges_C_per_s, ...
    min_state_exposure_s,min_plates_per_genotype_dTdt_bin, ...
    reorientation_cluster_gap_s,reorientation_cluster_mode, ...
    coldward_halfwidth_deg,warmward_threshold_deg, ...
    min_run_seconds,run_speed_cutoff_mm_s, ...
    top_directional_runs_to_plot,coldward_run_color,warmward_run_color, ...
    output_dir,"unc17_17_19C", ...
    "Figure6F_reorientation_vs_preceding_dTdt_WT_vs_unc17", ...
    "Figure6GH_longest_runs_WT_vs_unc17", ...
    "Supplement6FG_true_exit_rates_WT_vs_unc17", ...
    "Supplement6H_transition_matrices_WT_vs_unc17");

fprintf('\n============================================================\n');
fprintf('FIGURE 6 / SUPPLEMENT 6 ANALYSES COMPLETE\n');
fprintf('Outputs saved to:\n  %s\n',output_dir);
fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTIONS ADDED FOR THIS GITHUB-READY REDUCED SCRIPT
% ========================================================================

function thermalTrackT = local_extract_signed_thermal_velocity_tracks_only( ...
    coordsByStrain,plateIDdata,strainNames, ...
    xEdges,tempBinLabels,tempRange,pixRange,samplingRate)
% Minimal extractor for the signed thermal-migration metric used in Fig. 6A-B.
%
% This is the thermal-velocity portion of the larger original kinematic
% extractor, retained without the unrelated speed/pause metrics.
%
% Every consecutive valid processed coordinate pair contributes:
%
%   signedThermalVelocity = dT/dt
%
% X increases toward warm, so:
%   negative = coldward/cooling
%   positive = warmward/warming
%
% Each interval is assigned to a temperature bin by midpoint X, then interval
% velocities are averaged within each track x temperature bin.

    thermalTrackT = table();

    CperPx = diff(tempRange)/diff(pixRange);
    nBins = numel(tempBinLabels);

    for s = 1:numel(strainNames)
        coords = coordsByStrain{s};

        for w = 1:numel(coords)
            sumThermal = zeros(nBins,1);
            nThermal = zeros(nBins,1);

            if ~isfield(coords,'x_smooth') || isempty(coords(w).x_smooth)
                continue
            end

            xcell = coords(w).x_smooth;
            ycell = coords(w).y_smooth;
            if ~iscell(xcell), xcell = {xcell}; end
            if ~iscell(ycell), ycell = {ycell}; end

            for seg = 1:min(numel(xcell),numel(ycell))
                x = xcell{seg};
                y = ycell{seg};

                if isempty(x) || isempty(y)
                    continue
                end

                x = x(:);
                y = y(:);
                N = min(numel(x),numel(y));
                x = x(1:N);
                y = y(1:N);

                if N < 2
                    continue
                end

                x1 = x(1:end-1);
                x2 = x(2:end);
                y1 = y(1:end-1);
                y2 = y(2:end);

                good = isfinite(x1) & isfinite(x2) & ...
                       isfinite(y1) & isfinite(y2);

                if ~any(good)
                    continue
                end

                dx = x2(good)-x1(good);
                xMid = (x1(good)+x2(good))/2;

                signedThermalVel = ...
                    (dx.*CperPx).*samplingRate.*60;  % C/min

                binIdx = discretize(xMid,xEdges);

                for b = 1:nBins
                    kb = binIdx==b & isfinite(signedThermalVel);
                    if ~any(kb)
                        continue
                    end

                    sumThermal(b) = sumThermal(b) + ...
                        sum(signedThermalVel(kb),'omitnan');
                    nThermal(b) = nThermal(b) + ...
                        sum(isfinite(signedThermalVel(kb)));
                end
            end

            plateID = string(plateIDdata{s}(w));

            for b = 1:nBins
                if nThermal(b) < 1
                    continue
                end

                toks = regexp(tempBinLabels{b}, ...
                    '([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)', ...
                    'tokens','once');

                tempCenter = NaN;
                if ~isempty(toks)
                    tempCenter = mean(str2double(toks));
                end

                row = table( ...
                    string(strainNames{s}),plateID,w,b, ...
                    string(tempBinLabels{b}),tempCenter, ...
                    sumThermal(b)/nThermal(b), ...
                    'VariableNames',{ ...
                    'Strain','PlateID','TrackID','XBinIndex','XBin', ...
                    'TemperatureCenter_C', ...
                    'TrackMeanSignedThermalVelocity_C_per_min'});

                thermalTrackT = local_append_table(thermalTrackT,row);
            end
        end
    end
end


function R = local_run_paper_focused_roi( ...
    coordsByStrain,plateIDdata,plateOrder, ...
    strainNames,strainColors,compareIdx,roiC, ...
    samplingRate,pixRange,tempRange, ...
    dTdtWindow_s,headingWindow_s,dTdtEdges, ...
    minExposure_s,minPlatesPerDtdtBin, ...
    reversalGap_s,clusterMode,coldHalfWidth,warmThreshold, ...
    minRun_s,runSpeedCutoff,topDirectionalRunsN, ...
    coldRunColor,warmRunColor, ...
    outputDir,prefix, ...
    dTdtFigureName,runFigureName,exitFigureName,transitionFigureName)
% Reduced focused-ROI analysis containing only manuscript outcomes:
%   1) reorientation rate vs preceding 2-s dT/dt
%   2) longest coldward / warmward run visualization
%   3) true coldward-exit rate
%   4) true warmward-exit rate
%   5) incoming -> outgoing reorientation transition matrix

    compareNames = string(strainNames(compareIdx));

    [eventT,exposureT] = local_focus_build_reorientation_data( ...
        coordsByStrain,plateIDdata,plateOrder,strainNames,samplingRate, ...
        pixRange,tempRange,roiC,dTdtWindow_s,headingWindow_s,dTdtEdges, ...
        reversalGap_s,clusterMode,coldHalfWidth,warmThreshold);

    % ------------------------------------------------------------------
    % Figure 6C / 6F: clustered reorientation rate vs preceding 2-s dT/dt.
    % ------------------------------------------------------------------
    preRateT = local_focus_reorientation_rate_by_dTdt( ...
        eventT,exposureT,plateOrder,strainNames,dTdtEdges,minExposure_s,"Pre");

    preFit = local_focus_fit_rate_tuning_glmm( ...
        preRateT,compareNames,dTdtEdges,minPlatesPerDtdtBin);

    fig = local_plot_paper_reorientation_rate_vs_preceding_dTdt( ...
        preRateT,preFit,compareNames,strainNames,strainColors,roiC,dTdtEdges);

    local_save_clean_figure(fig,fullfile(outputDir,dTdtFigureName+".pdf"));

    % ------------------------------------------------------------------
    % Figure 6D-E / 6G-H: longest complete coldward / warmward runs.
    % ------------------------------------------------------------------
    runT = local_focus_build_clustered_complete_runs( ...
        coordsByStrain,plateIDdata,strainNames,samplingRate,pixRange,tempRange, ...
        roiC,reversalGap_s,clusterMode,coldHalfWidth,warmThreshold, ...
        minRun_s,runSpeedCutoff);

    [runFig,selectedRunsT] = local_focus_plot_top_directional_runs( ...
        runT,coordsByStrain,compareNames,strainNames,samplingRate, ...
        pixRange,tempRange,roiC,topDirectionalRunsN, ...
        coldRunColor,warmRunColor);

    local_save_clean_figure(runFig,fullfile(outputDir,runFigureName+".pdf"));

    % ------------------------------------------------------------------
    % Supplement 6C-D / 6F-G: TRUE directional exit rates.
    %
    % Coldward exit:
    %   #(coldward incoming -> transverse/warmward outgoing)
    %   / coldward animal-minutes
    %
    % Warmward exit:
    %   #(warmward incoming -> transverse/coldward outgoing)
    %   / warmward animal-minutes
    % ------------------------------------------------------------------
    coldExitT = local_focus_true_direction_exit_rate_table( ...
        eventT,exposureT,plateOrder,strainNames,minExposure_s,"Coldward");

    warmExitT = local_focus_true_direction_exit_rate_table( ...
        eventT,exposureT,plateOrder,strainNames,minExposure_s,"Warmward");

    coldExitFit = local_focus_fit_rate_glmm(coldExitT,compareNames);
    warmExitFit = local_focus_fit_rate_glmm(warmExitT,compareNames);

    coldExitPlate = local_focus_subset_plate_rate(coldExitT,compareNames);
    warmExitPlate = local_focus_subset_plate_rate(warmExitT,compareNames);

    exitFig = local_plot_paper_true_exit_rates( ...
        coldExitPlate,coldExitFit,warmExitPlate,warmExitFit, ...
        compareNames,strainNames,strainColors,roiC);

    local_save_clean_figure(exitFig,fullfile(outputDir,exitFigureName+".pdf"));

    % ------------------------------------------------------------------
    % Supplement 6E / 6H: transition matrices.
    % ------------------------------------------------------------------
    [transitionFig,transitionPlateT,transitionSummaryT] = ...
        local_focus_plot_transition_matrices( ...
        eventT,compareNames,strainNames,strainColors,roiC);

    local_save_clean_figure( ...
        transitionFig,fullfile(outputDir,transitionFigureName+".pdf"));

    % ------------------------------------------------------------------
    % Source-data exports.
    % ------------------------------------------------------------------
    pairMask = ismember(string(eventT.Strain),compareNames);
    eventPairT = eventT(pairMask,:);

    exposurePairT = exposureT( ...
        ismember(string(exposureT.Strain),compareNames),:);

    preRatePairT = preRateT( ...
        ismember(string(preRateT.Strain),compareNames),:);

    runPairT = runT(ismember(string(runT.Strain),compareNames),:);

    writetable(eventPairT,fullfile(outputDir,prefix+"_reorientation_events.csv"));
    writetable(exposurePairT,fullfile(outputDir,prefix+"_exposure.csv"));
    writetable(preRatePairT,fullfile(outputDir,prefix+"_preceding_dTdt_rates.csv"));
    writetable(preFit.BinContrasts, ...
        fullfile(outputDir,prefix+"_preceding_dTdt_bin_statistics.csv"));
    writetable(runPairT,fullfile(outputDir,prefix+"_complete_runs.csv"));
    writetable(selectedRunsT,fullfile(outputDir,prefix+"_selected_longest_runs.csv"));
    writetable(coldExitT(ismember(string(coldExitT.Strain),compareNames),:), ...
        fullfile(outputDir,prefix+"_true_coldward_exit_rates.csv"));
    writetable(warmExitT(ismember(string(warmExitT.Strain),compareNames),:), ...
        fullfile(outputDir,prefix+"_true_warmward_exit_rates.csv"));
    writetable(transitionPlateT, ...
        fullfile(outputDir,prefix+"_transition_probabilities_by_plate.csv"));
    writetable(transitionSummaryT, ...
        fullfile(outputDir,prefix+"_transition_matrix_summary.csv"));

    % Compact statistical summary.
    Analysis = [ ...
        "Genotype x preceding dT/dt"; ...
        "True coldward-exit rate"; ...
        "True warmward-exit rate"];

    Estimate = [NaN;coldExitFit.Estimate;warmExitFit.Estimate];
    CI95Low = [NaN;coldExitFit.CI95Low;warmExitFit.CI95Low];
    CI95High = [NaN;coldExitFit.CI95High;warmExitFit.CI95High];
    P = [preFit.InteractionP;coldExitFit.P;warmExitFit.P];

    Effect = [ ...
        "Genotype x dT/dt-bin interaction"; ...
        "Rate ratio mutant / WT"; ...
        "Rate ratio mutant / WT"];

    summaryT = table(Analysis,Effect,Estimate,CI95Low,CI95High,P);
    writetable(summaryT,fullfile(outputDir,prefix+"_statistics_summary.csv"));

    R = struct();
    R.ROI_C = roiC;
    R.CompareNames = compareNames;
    R.EventTable = eventPairT;
    R.ExposureTable = exposurePairT;
    R.PrecedingDtdtRateTable = preRatePairT;
    R.PrecedingDtdtFit = preFit;
    R.RunTable = runPairT;
    R.SelectedLongestRuns = selectedRunsT;
    R.ColdwardExitTable = coldExitT;
    R.ColdwardExitFit = coldExitFit;
    R.WarmwardExitTable = warmExitT;
    R.WarmwardExitFit = warmExitFit;
    R.TransitionMatrixPlateTable = transitionPlateT;
    R.TransitionMatrixSummaryTable = transitionSummaryT;
    R.Summary = summaryT;
end


function fig = local_plot_paper_thermal_velocity_pair( ...
    plateT,pairResult,strainNames,strainColors,compareIdx,minPlateN)
% Standalone curve corresponding to Figure 6A/B.
% One gray/colored point = one independent plate.
% Line/marker = mean across plates; error bars = SEM.
% Bin P values are Holm-corrected across eligible temperature bins.

    compareNames = string(strainNames(compareIdx));
    pairColors = strainColors(compareIdx,:);

    B = local_force_table(pairResult.BinContrasts);

    if ismember('TestStatus',B.Properties.VariableNames)
        B = B(string(B.TestStatus)=="tested",:);
    end

    B = B(isfinite(B.TemperatureCenter_C),:);
    B = sortrows(B,'TemperatureCenter_C','ascend');

    testedCenters = unique(B.TemperatureCenter_C,'sorted');

    fig = figure('Color','w','Position',[120 120 690 520]);
    ax = axes(fig);
    hold(ax,'on');

    dodge = [-0.06 0.06];

    yMaxData = -Inf;

    for gg = 1:2
        xSummary = [];
        ySummary = [];
        semSummary = [];

        for jj = 1:numel(testedCenters)
            c = testedCenters(jj);

            rows = abs(plateT.TemperatureCenter_C-c)<1e-9 & ...
                string(plateT.Strain)==compareNames(gg);

            vals = plateT.PlateMean(rows);
            vals = vals(isfinite(vals));

            if numel(vals) < minPlateN
                continue
            end

            x0 = c + dodge(gg);

            if numel(vals)==1
                jitter = 0;
            else
                jitter = linspace(-0.025,0.025,numel(vals))';
            end

            scatter(ax,x0+jitter,vals,26, ...
                'MarkerFaceColor',pairColors(gg,:), ...
                'MarkerEdgeColor','none', ...
                'MarkerFaceAlpha',0.25, ...
                'HandleVisibility','off');

            xSummary(end+1,1) = x0; %#ok<AGROW>
            ySummary(end+1,1) = mean(vals,'omitnan'); %#ok<AGROW>
            semSummary(end+1,1) = std(vals,0,'omitnan')/sqrt(numel(vals)); %#ok<AGROW>

            yMaxData = max(yMaxData,max(vals,[],'omitnan'));
        end

        errorbar(ax,xSummary,ySummary,semSummary,'-o', ...
            'Color',pairColors(gg,:), ...
            'MarkerFaceColor',pairColors(gg,:), ...
            'MarkerEdgeColor',pairColors(gg,:), ...
            'LineWidth',1.8, ...
            'MarkerSize',5, ...
            'CapSize',5, ...
            'DisplayName',char(compareNames(gg)));
    end

    yline(ax,0,'--','Color',[0.45 0.45 0.45], ...
        'LineWidth',0.9,'HandleVisibility','off');

    % Holm-corrected bin labels.
    yl = ylim(ax);
    yrange = diff(yl);
    if ~isfinite(yrange) || yrange<=0
        yrange = 1;
    end

    yText = yl(2)-0.05*yrange;

    for rr = 1:height(B)
        if ismember('p_Holm',B.Properties.VariableNames)
            p = B.p_Holm(rr);
        else
            p = B.p_raw(rr);
        end

        if ~isfinite(p)
            label = "";
        elseif p < 1e-4
            label = "p<0.0001";
        elseif p < 0.05
            label = "p="+string(sprintf('%.4f',p));
        else
            label = "ns";
        end

        text(ax,B.TemperatureCenter_C(rr),yText,label, ...
            'HorizontalAlignment','center','FontSize',9);
    end

    [~,~,pInt] = local_extract_lrt_result(pairResult.Interaction);

    text(ax,0.03,0.06,sprintf('Genotype x T\np = %s',local_pvalue_text(pInt)), ...
        'Units','normalized','VerticalAlignment','bottom');

    xlabel(ax,'Temperature (°C)');
    ylabel(ax,'Thermal velocity (°C/min)');
    legend(ax,'Location','best','Box','off');

    box(ax,'off');
    set(ax,'TickDir','out');
    ax.YGrid='on';
    ax.GridAlpha=0.15;

    title(ax,sprintf('%s vs %s',compareNames(1),compareNames(2)), ...
        'Interpreter','none','FontWeight','normal');
end


function fig = local_plot_paper_reorientation_rate_vs_preceding_dTdt( ...
    T,fit,compareNames,strainNames,strainColors,roiC,dTdtEdges)
% Standalone curve corresponding to Figure 6C/F.
%
% Each plate contributes one rate for each eligible preceding-dT/dt bin:
% clustered onset count / qualifying animal-minutes.
%
% Only dT/dt bins represented by >=4 independent plates in BOTH genotypes
% after the 10-s exposure filter are shown.

    fig = figure('Color','w','Position',[130 130 610 500]);
    ax = axes(fig);
    hold(ax,'on');

    validBins = [];
    if isstruct(fit) && isfield(fit,'AnalysisData') && ...
            istable(fit.AnalysisData) && ~isempty(fit.AnalysisData)
        validBins = sort(unique(fit.AnalysisData.BinIndex));
    end

    if isempty(validBins)
        text(ax,0.5,0.5,'No dT/dt bins meet analysis criteria', ...
            'Units','normalized','HorizontalAlignment','center');
        return
    end

    TT = T(ismember(T.BinIndex,validBins),:);

    for gg = 1:2
        sIdx = find(string(strainNames)==string(compareNames(gg)),1);

        centers = [];
        mu = [];
        se = [];

        for bb = validBins(:)'
            q = string(TT.Strain)==string(compareNames(gg)) & ...
                TT.BinIndex==bb & TT.Included;

            vals = TT.Rate(q);
            vals = vals(isfinite(vals));

            if isempty(vals)
                continue
            end

            centers(end+1,1) = mean(dTdtEdges(bb:bb+1)); %#ok<AGROW>
            mu(end+1,1) = mean(vals,'omitnan'); %#ok<AGROW>

            if numel(vals)>=2
                se(end+1,1) = std(vals,0,'omitnan')/sqrt(numel(vals)); %#ok<AGROW>
            else
                se(end+1,1) = NaN; %#ok<AGROW>
            end
        end

        errorbar(ax,centers,mu,se,'-o', ...
            'Color',strainColors(sIdx,:), ...
            'MarkerFaceColor',strainColors(sIdx,:), ...
            'MarkerEdgeColor',strainColors(sIdx,:), ...
            'LineWidth',1.8,'MarkerSize',5,'CapSize',5, ...
            'DisplayName',char(compareNames(gg)));
    end

    xline(ax,0,':','Color',[0.55 0.55 0.55], ...
        'HandleVisibility','off');

    xlabel(ax,'Preceding 2-s mean dT/dt (C/s)');
    ylabel(ax,'Reorientation rate');
    xlim(ax,[dTdtEdges(1) dTdtEdges(end)]);
    xticks(ax,dTdtEdges);

    text(ax,0.03,0.94,sprintf('Genotype x dT/dt\np = %s', ...
        local_pvalue_text(fit.InteractionP)), ...
        'Units','normalized','VerticalAlignment','top');

    title(ax,sprintf('%.0f-%.0f °C',roiC(1),roiC(2)), ...
        'FontWeight','normal');

    legend(ax,'Location','best','Box','off');
    box(ax,'off');
    set(ax,'TickDir','out');
    ax.YGrid='on';
    ax.GridAlpha=0.15;
end


function fig = local_plot_paper_true_exit_rates( ...
    coldPlate,coldFit,warmPlate,warmFit, ...
    compareNames,strainNames,strainColors,roiC)
% Two-panel source figure corresponding to Supplement 6 C-D or F-G.

    fig = figure('Color','w','Position',[100 120 900 430]);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    local_focus_plate_dot_panel( ...
        nexttile(tl,1),coldPlate,compareNames,strainNames,strainColors, ...
        'Coldward exit rate','Exits / coldward animal-min',coldFit,NaN);

    local_focus_plate_dot_panel( ...
        nexttile(tl,2),warmPlate,compareNames,strainNames,strainColors, ...
        'Warmward exit rate','Exits / warmward animal-min',warmFit,NaN);

    title(tl,sprintf('%s vs %s | %.0f-%.0f °C', ...
        compareNames(1),compareNames(2),roiC(1),roiC(2)), ...
        'Interpreter','none','FontWeight','normal');
end


function s = local_pvalue_text(p)
    if ~isfinite(p)
        s = 'NA';
    elseif p < 1e-4
        s = '<0.0001';
    elseif p < 0.001
        s = sprintf('%.2g',p);
    elseif p < 0.01
        s = sprintf('%.4f',p);
    else
        s = sprintf('%.4f',p);
    end
end


function T = local_sort_temperature_table(T)
% Force deterministic numeric temperature-bin order.
%
% Never rely on first-observed row order: individual tracks can visit only a
% subset of bins, so a track-major table can first encounter bins in an order
% such as 15-16, 17-18, 16-17. That is valid data but not valid plotting order.

    if isempty(T) || ~istable(T)
        return
    end

    vars=T.Properties.VariableNames;

    % Validate XBinIndex -> TemperatureCenter_C mapping when both exist.
    if ismember('XBinIndex',vars) && ismember('TemperatureCenter_C',vars)
        q=isfinite(T.XBinIndex) & isfinite(T.TemperatureCenter_C);

        if any(q)
            bins=unique(T.XBinIndex(q));
            centers=nan(numel(bins),1);

            for ii=1:numel(bins)
                vv=T.TemperatureCenter_C(q & T.XBinIndex==bins(ii));
                centers(ii)=mean(vv,'omitnan');

                if numel(vv)>1 && (max(vv)-min(vv))>1e-9
                    error('Inconsistent TemperatureCenter_C within XBinIndex %g.',bins(ii));
                end
            end

            [binsSorted,ord]=sort(bins); %#ok<ASGLU>
            centersSorted=centers(ord);

            if numel(centersSorted)>1 && any(diff(centersSorted)<=0)
                error(['TemperatureCenter_C does not increase with XBinIndex. ' ...
                       'The bin-to-temperature mapping is inconsistent.']);
            end
        end
    end

    sortVars={};
    sortDirs={};

    if ismember('Strain',vars)
        sortVars{end+1}='Strain'; %#ok<AGROW>
        sortDirs{end+1}='ascend'; %#ok<AGROW>
    end
    if ismember('PlateID',vars)
        sortVars{end+1}='PlateID'; %#ok<AGROW>
        sortDirs{end+1}='ascend'; %#ok<AGROW>
    end
    if ismember('TrackID',vars)
        sortVars{end+1}='TrackID'; %#ok<AGROW>
        sortDirs{end+1}='ascend'; %#ok<AGROW>
    end

    if ismember('XBinIndex',vars)
        sortVars{end+1}='XBinIndex'; %#ok<AGROW>
        sortDirs{end+1}='ascend'; %#ok<AGROW>
    elseif ismember('TemperatureCenter_C',vars)
        sortVars{end+1}='TemperatureCenter_C'; %#ok<AGROW>
        sortDirs{end+1}='ascend'; %#ok<AGROW>
    end

    if ~isempty(sortVars)
        T=sortrows(T,sortVars,sortDirs);
    end
end


function T = local_append_table(T,row)
% Safely append one or more rows to a table that may initially be empty.
% Several analyses in this script build long-format tables incrementally.
%
% Using this helper avoids MATLAB errors when concatenating into table().
    if isempty(T) || width(T)==0
        T = row;
    else
        T = [T; row]; %#ok<AGROW>
    end
end


function plateCurveT = local_track_table_to_plate_curve_summary(trackT,metricVar)
% Mean of track-level metric values within each independent plate x bin.

    plateCurveT=table();
    if isempty(trackT),return;end

    keys=unique(trackT(:,{'Strain','PlateID','XBinIndex','XBin','TemperatureCenter_C'}), ...
        'rows','stable');

    for k=1:height(keys)
        rows=string(trackT.Strain)==string(keys.Strain(k)) & ...
            string(trackT.PlateID)==string(keys.PlateID(k)) & ...
            trackT.XBinIndex==keys.XBinIndex(k);

        vals=trackT.(metricVar)(rows);
        vals=vals(isfinite(vals));

        if isempty(vals)
            plateMean=NaN;
            nTracks=0;
        else
            plateMean=mean(vals,'omitnan');
            nTracks=numel(vals);
        end

        row=table(string(keys.Strain(k)),string(keys.PlateID(k)), ...
            keys.XBinIndex(k),string(keys.XBin(k)),keys.TemperatureCenter_C(k), ...
            plateMean,nTracks, ...
            'VariableNames',{'Strain','PlateID','XBinIndex','XBin', ...
            'TemperatureCenter_C','PlateMean','NTracks'});
        plateCurveT=local_append_table(plateCurveT,row);
    end

    % First-observed row order is not guaranteed to match temperature order.
    plateCurveT=local_sort_temperature_table(plateCurveT);
end


function rowT = local_lme_residual_diagnostic_row(lme,scaleName)
% Compact residual diagnostics for a fitted LME.
% These are diagnostic summaries, not independent hypothesis tests used to
% choose significance. Smaller absolute skew/kurtosis and weaker association
% between |residual| and fitted value indicate improved model behavior.

    try
        r=residuals(lme,'Conditional',true);
    catch
        r=residuals(lme);
    end
    try
        f=fitted(lme,'Conditional',true);
    catch
        f=fitted(lme);
    end

    r=r(:); f=f(:);
    keep=isfinite(r) & isfinite(f);
    r=r(keep); f=f(keep);

    N=numel(r);
    Scale=string(scaleName);

    if N>=3
        ResidualSkewness=skewness(r,0);
        AbsResidualSkewness=abs(ResidualSkewness);
        ExcessKurtosis=kurtosis(r,0)-3;
        AbsExcessKurtosis=abs(ExcessKurtosis);
    else
        ResidualSkewness=NaN;
        AbsResidualSkewness=NaN;
        ExcessKurtosis=NaN;
        AbsExcessKurtosis=NaN;
    end

    % Jarque-Bera residual normality diagnostic.
    %
    % Do NOT call MATLAB jbtest here. jbtest uses tabulated p-values and,
    % for very non-normal residuals, emits:
    %   "P is less than the smallest tabulated value, returning 0.001."
    %
    % For this diagnostic comparison we calculate the standard asymptotic
    % Jarque-Bera statistic directly:
    %
    %   JB = N/6 * [skewness^2 + excess_kurtosis^2/4]
    %
    % Under the null, JB is asymptotically chi-square with 2 df. For 2 df,
    % the upper-tail probability is exactly exp(-JB/2), so this gives a
    % continuous diagnostic p-value without the tabulation warning.
    JarqueBera_stat=NaN;
    JarqueBera_p=NaN;
    JarqueBera_h=NaN;

    if N>=4 && isfinite(ResidualSkewness) && isfinite(ExcessKurtosis)
        JarqueBera_stat = (N/6) * ...
            (ResidualSkewness^2 + (ExcessKurtosis^2)/4);

        JarqueBera_p = exp(-JarqueBera_stat/2);
        JarqueBera_h = JarqueBera_p < 0.05;
    end

    SpearmanRho_AbsResidual_vs_Fitted=NaN;
    SpearmanP_AbsResidual_vs_Fitted=NaN;
    if N>=3 && numel(unique(f))>1
        try
            [rho,pv]=corr(f,abs(r),'Type','Spearman','Rows','complete');
            SpearmanRho_AbsResidual_vs_Fitted=rho;
            SpearmanP_AbsResidual_vs_Fitted=pv;
        catch
            C=corrcoef(f,abs(r),'Rows','complete');
            if all(size(C)==[2 2])
                SpearmanRho_AbsResidual_vs_Fitted=C(1,2);
            end
        end
    end
    AbsSpearmanRho_AbsResidual_vs_Fitted= ...
        abs(SpearmanRho_AbsResidual_vs_Fitted);

    ResidualSD=std(r,0,'omitnan');

    rowT=table(Scale,N,ResidualSD,ResidualSkewness,AbsResidualSkewness, ...
        ExcessKurtosis,AbsExcessKurtosis,JarqueBera_stat,JarqueBera_h,JarqueBera_p, ...
        SpearmanRho_AbsResidual_vs_Fitted, ...
        AbsSpearmanRho_AbsResidual_vs_Fitted, ...
        SpearmanP_AbsResidual_vs_Fitted);
end


function [lmeFull,interactionCompare,overallGenotypeT,binContrastT,modelInfoT, ...
    fixedEffectsT,diagnosticT,analysisT] = ...
    local_fit_signed_track_metric_lme( ...
    trackT,strainNames,compareIdx,tempLabels,minPlatesPerGroup, ...
    multipleComparisonMethod,metricName)
% Hierarchical LME for signed TRACK-level measurements such as DeltaAlignment.
%
% Response ~ Genotype*TemperatureBin + (1|PlateUID) + (1|TrackUID)
%
% Track observations are retained, but plate remains the biological replicate.

    compareNames=string(strainNames(compareIdx));
    tempLabelsS=string(tempLabels(:));
    multipleComparisonMethod=lower(string(multipleComparisonMethod));

    if multipleComparisonMethod~="holm"
        error('Only Holm multiple-comparison correction is currently implemented.');
    end

    metricVars=setdiff(trackT.Properties.VariableNames, ...
        {'Strain','PlateID','TrackID','XBinIndex','XBin','TemperatureCenter_C'});
    if numel(metricVars)~=1
        error('Could not uniquely identify signed track metric column.');
    end
    metricVar=metricVars{1};

    keep=ismember(string(trackT.Strain),compareNames) & ...
         isfinite(trackT.(metricVar));
    T=trackT(keep,:);

    if isempty(T)
        error('No finite track-level observations remain for %s.',metricName);
    end

    T.PlateUID=categorical(string(T.Strain)+"__plate__"+string(T.PlateID));
    T.TrackUID=categorical(string(T.Strain)+"__plate__"+string(T.PlateID)+ ...
        "__track__"+string(T.TrackID));

    nBins=numel(tempLabelsS);
    nControlPlates=zeros(nBins,1);
    nMutantPlates=zeros(nBins,1);
    nControlTracks=zeros(nBins,1);
    nMutantTracks=zeros(nBins,1);
    modelBinMask=false(nBins,1);

    % GLOBAL MODEL BIN ELIGIBILITY:
    % Require the same minimum independent-plate support as the bin contrasts.
    for b=1:nBins
        inBin=string(T.XBin)==tempLabelsS(b);
        kc=inBin & string(T.Strain)==compareNames(1);
        km=inBin & string(T.Strain)==compareNames(2);

        nControlPlates(b)=numel(unique(string(T.PlateUID(kc))));
        nMutantPlates(b)=numel(unique(string(T.PlateUID(km))));
        nControlTracks(b)=sum(kc);
        nMutantTracks(b)=sum(km);
        modelBinMask(b)=nControlPlates(b)>=minPlatesPerGroup && ...
            nMutantPlates(b)>=minPlatesPerGroup;
    end

    modelBinLabels=tempLabelsS(modelBinMask);
    if numel(modelBinLabels)<2
        error(['%s hierarchical LME requires at least two temperature bins with ' ...
            '>= %d independent plates in BOTH genotypes.'],metricName,minPlatesPerGroup);
    end

    T=T(ismember(string(T.XBin),modelBinLabels),:);
    T=local_sort_temperature_table(T);
    T.Genotype=categorical(string(T.Strain),compareNames,compareNames);
    T.TemperatureBin=categorical(string(T.XBin),modelBinLabels,modelBinLabels);
    T.Response=T.(metricVar);

    analysisT=T(:,{'Strain','PlateID','PlateUID','TrackID','TrackUID', ...
        'XBinIndex','XBin','TemperatureCenter_C','Response'});

    fullFormula=[ ...
        'Response ~ Genotype*TemperatureBin + ' ...
        '(1|PlateUID) + (1|TrackUID)'];

    reducedFormula=[ ...
        'Response ~ Genotype + TemperatureBin + ' ...
        '(1|PlateUID) + (1|TrackUID)'];

    temperatureOnlyFormula=[ ...
        'Response ~ TemperatureBin + ' ...
        '(1|PlateUID) + (1|TrackUID)'];

    lmeFull=fitlme(T,fullFormula,'FitMethod','ML','DummyVarCoding','reference');
    lmeReduced=fitlme(T,reducedFormula,'FitMethod','ML','DummyVarCoding','reference');
    lmeTemperatureOnly=fitlme(T,temperatureOnlyFormula, ...
        'FitMethod','ML','DummyVarCoding','reference');

    % Global Genotype x Temperature interaction.
    interactionCompare=local_force_table( ...
        compare(lmeReduced,lmeFull,'CheckNesting',true));

    % Overall genotype effect across temperature from the additive model.
    overallLRT=local_force_table( ...
        compare(lmeTemperatureOnly,lmeReduced,'CheckNesting',true));

    betaAdd=fixedEffects(lmeReduced);
    VAdd=lmeReduced.CoefficientCovariance;
    XAdd=designMatrix(lmeReduced,'Fixed');

    firstBin=modelBinLabels(1);
    iControlAdd=find(string(T.XBin)==firstBin & ...
        string(T.Strain)==compareNames(1),1,'first');
    iMutantAdd=find(string(T.XBin)==firstBin & ...
        string(T.Strain)==compareNames(2),1,'first');

    HAdd=XAdd(iMutantAdd,:)-XAdd(iControlAdd,:);
    overallDiff=HAdd*betaAdd;
    overallSE=sqrt(max(real(HAdd*VAdd*HAdd'),0));

    overallDF=NaN;
    overallWaldP=NaN;
    try
        [overallWaldP,~,~,overallDF]=coefTest( ...
            lmeReduced,HAdd,0,'DFMethod','Satterthwaite');
    catch
        try
            [overallWaldP,~,~,overallDF]=coefTest(lmeReduced,HAdd,0);
        catch
            if isfinite(overallSE) && overallSE>0
                overallWaldP=erfc(abs(overallDiff/overallSE)/sqrt(2));
            end
        end
    end

    if isfinite(overallDF) && overallDF>0
        overallCrit=tinv(0.975,overallDF);
    else
        overallCrit=1.959963984540054;
    end

    [overallLRStat,overallDeltaDF,overallLRTP] = ...
        local_extract_lrt_result(overallLRT);

    ControlGenotypeOverall=compareNames(1);
    MutantGenotypeOverall=compareNames(2);
    ModelAssumption="common genotype effect across temperature (additive model)";
    N_control_plates_total=numel(unique(string(T.PlateUID(string(T.Strain)==compareNames(1)))));
    N_mutant_plates_total=numel(unique(string(T.PlateUID(string(T.Strain)==compareNames(2)))));
    N_control_tracks_total=numel(unique(string(T.TrackUID(string(T.Strain)==compareNames(1)))));
    N_mutant_tracks_total=numel(unique(string(T.TrackUID(string(T.Strain)==compareNames(2)))));
    Difference_mutant_minus_control=overallDiff;
    Difference_CI95_low=overallDiff-overallCrit*overallSE;
    Difference_CI95_high=overallDiff+overallCrit*overallSE;
    SE_Difference=overallSE;
    ContrastDF=overallDF;
    Wald_p=overallWaldP;
    LRT_LRStat=overallLRStat;
    LRT_deltaDF=overallDeltaDF;
    LRT_pValue=overallLRTP;
    Significant_LRT_0p05=isfinite(LRT_pValue) && LRT_pValue<0.05;

    overallGenotypeT=table( ...
        ControlGenotypeOverall,MutantGenotypeOverall,ModelAssumption, ...
        N_control_plates_total,N_mutant_plates_total, ...
        N_control_tracks_total,N_mutant_tracks_total, ...
        Difference_mutant_minus_control,Difference_CI95_low,Difference_CI95_high, ...
        SE_Difference,ContrastDF,Wald_p,LRT_LRStat,LRT_deltaDF,LRT_pValue, ...
        Significant_LRT_0p05);

    beta=fixedEffects(lmeFull);
    fixedEffectsT=local_force_table(lmeFull.Coefficients);
    V=lmeFull.CoefficientCovariance;
    Xfixed=designMatrix(lmeFull,'Fixed');

    diagnosticT=local_lme_residual_diagnostic_row(lmeFull,"untransformed_track_metric");

    TemperatureBin=tempLabelsS;
    TemperatureCenter_C=nan(nBins,1);
    N_control_plates=nControlPlates;
    N_mutant_plates=nMutantPlates;
    N_control_tracks=nControlTracks;
    N_mutant_tracks=nMutantTracks;
    ControlGenotype=repmat(compareNames(1),nBins,1);
    MutantGenotype=repmat(compareNames(2),nBins,1);

    ModelMean_control=nan(nBins,1);
    ModelMean_mutant=nan(nBins,1);
    Difference_mutant_minus_control=nan(nBins,1);
    Difference_CI95_low=nan(nBins,1);
    Difference_CI95_high=nan(nBins,1);
    SE_Difference=nan(nBins,1);
    ContrastDF=nan(nBins,1);
    p_raw=nan(nBins,1);
    p_Holm=nan(nBins,1);
    Significant_Holm_0p05=false(nBins,1);
    TestStatus=strings(nBins,1);

    for b=1:nBins
        toks=regexp(char(tempLabelsS(b)), ...
            '([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)','tokens','once');
        if ~isempty(toks)
            TemperatureCenter_C(b)=mean(str2double(toks));
        end

        if N_control_plates(b)<minPlatesPerGroup || ...
                N_mutant_plates(b)<minPlatesPerGroup
            TestStatus(b)="not tested: fewer than minimum independent plates in >=1 genotype";
            continue
        end

        if ~ismember(tempLabelsS(b),modelBinLabels)
            TestStatus(b)="not estimable in global hierarchical LME";
            continue
        end

        iControl=find(string(T.XBin)==tempLabelsS(b) & ...
            string(T.Strain)==compareNames(1),1,'first');
        iMutant=find(string(T.XBin)==tempLabelsS(b) & ...
            string(T.Strain)==compareNames(2),1,'first');

        if isempty(iControl) || isempty(iMutant)
            TestStatus(b)="not estimable in global hierarchical LME";
            continue
        end

        xControl=Xfixed(iControl,:);
        xMutant=Xfixed(iMutant,:);
        H=xMutant-xControl;

        diffVal=H*beta;
        varDiff=max(real(H*V*H'),0);
        seDiff=sqrt(varDiff);

        df2=NaN;
        try
            [pThis,~,~,df2]=coefTest(lmeFull,H,0,'DFMethod','Satterthwaite');
        catch
            try
                [pThis,~,~,df2]=coefTest(lmeFull,H,0);
            catch
                if isfinite(seDiff) && seDiff>0
                    pThis=erfc(abs(diffVal/seDiff)/sqrt(2));
                else
                    pThis=NaN;
                end
            end
        end

        if isfinite(df2) && df2>0
            crit=tinv(0.975,df2);
        else
            crit=1.959963984540054;
        end

        ModelMean_control(b)=xControl*beta;
        ModelMean_mutant(b)=xMutant*beta;
        Difference_mutant_minus_control(b)=diffVal;
        Difference_CI95_low(b)=diffVal-crit*seDiff;
        Difference_CI95_high(b)=diffVal+crit*seDiff;
        SE_Difference(b)=seDiff;
        ContrastDF(b)=df2;
        p_raw(b)=pThis;
        TestStatus(b)="tested";
    end

    tested=isfinite(p_raw) & TestStatus=="tested";
    if any(tested)
        p_Holm(tested)=local_holm_adjust(p_raw(tested));
        Significant_Holm_0p05(tested)=p_Holm(tested)<0.05;
    end

    binContrastT=table(TemperatureBin,TemperatureCenter_C, ...
        ControlGenotype,MutantGenotype, ...
        N_control_plates,N_mutant_plates,N_control_tracks,N_mutant_tracks, ...
        ModelMean_control,ModelMean_mutant,Difference_mutant_minus_control, ...
        Difference_CI95_low,Difference_CI95_high,SE_Difference,ContrastDF, ...
        p_raw,p_Holm,Significant_Holm_0p05,TestStatus);

    Model=string(metricName)+" | track-level hierarchical Gaussian LME";
    ObservationUnit="track x temperature bin";
    IndependentBiologicalReplicate="plate";
    RandomEffects="PlateUID + TrackUID";
    Transformation="none";
    N_observations=height(T);
    N_independent_plates=numel(unique(string(T.PlateUID)));
    N_unique_tracks=numel(unique(string(T.TrackUID)));
    N_temperature_bins_in_model=numel(modelBinLabels);
    MinimumIndependentPlatesPerGenotypeForBinContrast=minPlatesPerGroup;
    MultipleComparisonCorrection=multipleComparisonMethod;

    modelInfoT=table(Model,ObservationUnit,IndependentBiologicalReplicate, ...
        RandomEffects,Transformation,N_observations,N_independent_plates, ...
        N_unique_tracks,N_temperature_bins_in_model, ...
        MinimumIndependentPlatesPerGenotypeForBinContrast, ...
        MultipleComparisonCorrection);
end


function T = local_force_table(x)
% Convert common MATLAB tabular outputs to table.
    if istable(x)
        T=x;
        return
    end

    % fitlme/fitglme and compare can return legacy dataset arrays.
    try
        if isa(x,'dataset')
            T=dataset2table(x);
            return
        end
    catch
    end

    % Some releases expose objects convertible via struct.
    try
        sx=struct(x);
        if ~isempty(sx)
            T=struct2table(sx);
            return
        end
    catch
    end

    T=x;
end


function [LRStat,deltaDF,pValue] = local_extract_lrt_result(lrtT)
% Extract the likelihood-ratio statistic, delta DF, and p-value from the
% second (larger-model) row returned by compare().
    LRStat=NaN;
    deltaDF=NaN;
    pValue=NaN;

    lrtT=local_force_table(lrtT);
    if ~istable(lrtT) || height(lrtT)<1
        return
    end

    row=height(lrtT);

    if ismember('LRStat',lrtT.Properties.VariableNames)
        LRStat=lrtT.LRStat(row);
    end
    if ismember('deltaDF',lrtT.Properties.VariableNames)
        deltaDF=lrtT.deltaDF(row);
    end
    if ismember('pValue',lrtT.Properties.VariableNames)
        pValue=lrtT.pValue(row);
    end
end


function mu=local_circ_mean_deg(theta_deg)
    th=deg2rad(theta_deg(:));
    mu=rad2deg(atan2(mean(sin(th),'omitnan'),mean(cos(th),'omitnan')));
end


function padj=local_holm_adjust(p)
    padj=nan(size(p)); valid=find(isfinite(p));
    if isempty(valid); return; end
    pv=p(valid); [ps,ord]=sort(pv); m=numel(ps);
    a=zeros(m,1);
    for k=1:m; a(k)=min(1,(m-k+1)*ps(k)); end
    for k=2:m; a(k)=max(a(k),a(k-1)); end
    tmp=zeros(m,1); tmp(ord)=a; padj(valid)=tmp;
end


function [X,Y,plateID] = local_load_xy_and_plate(xFile,yFile,plateFile,xVar,yVar,plateCandidates)
% Load X, Y, and plate identity when they are stored in SEPARATE MAT files.
%
% Expected organization, for example:
%   X_N2.mat
%   Y_N2.mat
%   plate_N2.mat
%
% The plate MAT file must contain one plate-ID value per trajectory column.
% Variable-name matching is permissive/case-insensitive.

    % ------------------------------------------------------------
    % Check that all three files exist
    % ------------------------------------------------------------
    if ~isfile(xFile)
        error('X file not found: %s',xFile);
    end
    if ~isfile(yFile)
        error('Y file not found: %s',yFile);
    end
    if ~isfile(plateFile)
        error('Plate-ID file not found: %s',plateFile);
    end

    SX = load(xFile);
    SY = load(yFile);
    SP = load(plateFile);

    xVar = char(xVar);
    yVar = char(yVar);

    % ------------------------------------------------------------
    % Load X and Y
    % ------------------------------------------------------------
    if ~isfield(SX,xVar)
        error('Variable "%s" was not found in %s',xVar,xFile);
    end
    if ~isfield(SY,yVar)
        error('Variable "%s" was not found in %s',yVar,yFile);
    end

    X = SX.(xVar);
    Y = SY.(yVar);

    if size(X,2) ~= size(Y,2)
        error('%s and %s have different numbers of trajectory columns (%d vs %d).', ...
            xVar,yVar,size(X,2),size(Y,2));
    end

    nTracks = size(X,2);

    % ------------------------------------------------------------
    % Find plate vector INSIDE plateFile
    % ------------------------------------------------------------
    plateID = [];
    foundName = '';
    plateCandidates = cellstr(string(plateCandidates));
    fp = fieldnames(SP);

    % 1) Exact candidate-name match
    for k = 1:numel(plateCandidates)
        nm = plateCandidates{k};
        if isfield(SP,nm)
            candidate = SP.(nm);
            if local_is_valid_plate_vector(candidate,nTracks)
                plateID = candidate;
                foundName = nm;
                break
            end
        end
    end

    % 2) Case-insensitive candidate-name match
    if isempty(plateID)
        for k = 1:numel(plateCandidates)
            target = plateCandidates{k};
            ii = find(strcmpi(fp,target),1);
            if ~isempty(ii)
                nm = fp{ii};
                candidate = SP.(nm);
                if local_is_valid_plate_vector(candidate,nTracks)
                    plateID = candidate;
                    foundName = nm;
                    break
                end
            end
        end
    end

    % 3) Prefer any correctly-sized variable whose name contains "plate"
    if isempty(plateID)
        [plateID,foundName] = local_find_plate_vector_in_struct(SP,nTracks);
    end

    % 4) Last-resort: if plateFile contains exactly ONE correctly-sized vector,
    %    use it even if its variable name does not contain "plate".
    if isempty(plateID)
        validNames = {};
        for k = 1:numel(fp)
            candidate = SP.(fp{k});
            if local_is_valid_plate_vector(candidate,nTracks)
                validNames{end+1} = fp{k}; %#ok<AGROW>
            end
        end

        if numel(validNames) == 1
            foundName = validNames{1};
            plateID = SP.(foundName);
        elseif numel(validNames) > 1
            fprintf('\nVariables in %s with %d elements:\n',plateFile,nTracks);
            fprintf('  %s\n',validNames{:});
            error(['Multiple possible plate-ID vectors were found in %s. ' ...
                   'Rename the intended variable to include "plate" or add its name to plateCandidates.'], ...
                   plateFile);
        end
    end

    % ------------------------------------------------------------
    % Diagnostic error if plate ID still not found
    % ------------------------------------------------------------
    if isempty(plateID)
        fprintf('\nVariables found in %s:\n',plateFile);
        for k = 1:numel(fp)
            v = SP.(fp{k});
            sz = size(v);
            fprintf('  %-30s size = %s  class = %s\n', ...
                fp{k},mat2str(sz),class(v));
        end

        error(['No valid plate-ID vector could be identified in %s.\n' ...
               'X has %d trajectory columns, so the plate vector must contain exactly %d entries.'], ...
               plateFile,nTracks,nTracks);
    end

    % ------------------------------------------------------------
    % Normalize to one plate ID per trajectory column
    % ------------------------------------------------------------
    plateID = plateID(:);

    if numel(plateID) ~= nTracks
        error(['Plate-ID length mismatch for %s: X has %d trajectory columns, ' ...
               'but "%s" in %s has %d entries.'], ...
               xVar,nTracks,foundName,plateFile,numel(plateID));
    end

    fprintf(['Loaded %-12s | tracks = %d | plate variable = "%s" | ' ...
             'independent plates = %d\n'], ...
        xVar,nTracks,foundName,numel(unique(string(plateID))));
end


function tf = local_is_valid_plate_vector(v,nTracks)
% A valid plate variable has exactly one element per trajectory column.
    tf = (isnumeric(v) || islogical(v) || isstring(v) || iscellstr(v) || ...
          iscategorical(v) || ischar(v)) && numel(v)==nTracks;
end


function [plateID,foundName] = local_find_plate_vector_in_struct(S,nTracks)
% Find a plausible plate-ID variable automatically.
% Preference is given to names containing both "plate" and "id".

    plateID = [];
    foundName = '';
    fn = fieldnames(S);

    candidates = {};
    scores = [];

    for k = 1:numel(fn)
        nm = fn{k};
        low = lower(nm);

        if ~contains(low,'plate')
            continue
        end

        v = S.(nm);
        if ~local_is_valid_plate_vector(v,nTracks)
            continue
        end

        score = 1;
        if contains(low,'id')
            score = score + 2;
        end
        if strcmpi(nm,'plate_id') || strcmpi(nm,'plateid')
            score = score + 3;
        end

        candidates{end+1} = nm; %#ok<AGROW>
        scores(end+1) = score; %#ok<AGROW>
    end

    if isempty(candidates)
        return
    end

    [~,best] = max(scores);
    foundName = candidates{best};
    plateID = S.(foundName);
end


function clusters = local_cluster_reversal_calls(idx,gap_fr,clusterMode)
% Return [first_index last_index count].
%
% 'first' mode prevents chaining:
% a candidate joins only if it is within gap_fr of the FIRST call.
%
% 'previous' mode reproduces the old behavior:
% a candidate joins if it is within gap_fr of the PREVIOUS call.

    idx = idx(:);
    clusters = zeros(0,3);
    if isempty(idx), return; end

    firstK = idx(1);
    lastK = idx(1);
    count = 1;

    for k = 2:numel(idx)
        if clusterMode == "first"
            joinCluster = (idx(k)-firstK) <= gap_fr;
        else
            joinCluster = (idx(k)-lastK) <= gap_fr;
        end

        if joinCluster
            lastK = idx(k);
            count = count + 1;
        else
            clusters(end+1,:) = [firstK,lastK,count]; %#ok<AGROW>
            firstK = idx(k);
            lastK = idx(k);
            count = 1;
        end
    end

    clusters(end+1,:) = [firstK,lastK,count];
end


function heading = local_heading_segment_for_reversal(x,y)
% ============================================================
% DIAGNOSTIC VERSION
%
% Purpose:
% Test whether calculateRelativeHeadingandDistance ever fails during the
% actual reversal-heading analysis.
%
% Expected heading convention:
%   0 deg      = coldward
%   +/-90 deg  = isothermal
%   +/-180 deg = warmward
%
% IMPORTANT:
% There is intentionally NO atan2d fallback in this diagnostic version.
% If calculateRelativeHeadingandDistance fails, MATLAB will stop and show
% the original error.
% ============================================================

    x = x(:);
    y = y(:);

    try
        [heading,~] = calculateRelativeHeadingandDistance( ...
            x,y,[-1 0]);

    catch ME
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('ERROR IN calculateRelativeHeadingandDistance\n');
        fprintf('============================================================\n');
        fprintf('Number of X points: %d\n',numel(x));
        fprintf('Number of Y points: %d\n',numel(y));
        fprintf('Finite X points: %d\n',sum(isfinite(x)));
        fprintf('Finite Y points: %d\n',sum(isfinite(y)));
        fprintf('\nMATLAB error message:\n%s\n',ME.message);
        fprintf('============================================================\n\n');

        % Stop immediately and display the complete original MATLAB error.
        rethrow(ME)
    end

    % ------------------------------------------------------------
    % Make output length match trajectory length
    % ------------------------------------------------------------
    if numel(heading) == numel(x)-1

        % Repeat final heading value so there is one heading per coordinate.
        heading(end+1,1) = heading(end);

    elseif numel(heading) ~= numel(x)

        error(['calculateRelativeHeadingandDistance returned %d headings ' ...
               'for a trajectory containing %d coordinates.'], ...
               numel(heading),numel(x));
    end

    heading = heading(:);

    % ------------------------------------------------------------
    % Additional sanity check
    % ------------------------------------------------------------
    finiteHeading = heading(isfinite(heading));

    if ~isempty(finiteHeading)
        if any(finiteHeading < -180 | finiteHeading > 180)
            warning(['Unexpected heading values outside [-180,180] deg. ' ...
                     'Observed range = %.2f to %.2f deg.'], ...
                     min(finiteHeading),max(finiteHeading));
        end
    end
end


function [fig,plateProbT,summaryT] = local_focus_plot_transition_matrices( ...
    eventT,compareNames,strainNames,strainColors,roiC)
% ROI-SPECIFIC REORIENTATION TRANSITION MATRICES.
%
% Rows    = incoming direction before the clustered reorientation.
% Columns = outgoing direction after the clustered reorientation.
%
% The matrix answers:
%
%   P(outgoing direction | incoming direction, a reorientation occurred)
%
% IMPORTANT HIERARCHY:
%   1) Within each independent plate and incoming direction, count outgoing
%      Coldward / Transverse / Warmward choices.
%   2) Divide by that plate's total number of valid transitions from that
%      incoming direction so the three probabilities sum to 1.
%   3) Average those plate-level probabilities across plates.
%
% Therefore each displayed ROW sums to 1 (apart from rounding), while each
% plate receives equal weight. Plates with zero valid transitions for a given
% incoming direction have an undefined conditional distribution and do not
% contribute to that row.
%
% Only events with BOTH a valid IncomingDirection and OutgoingDirection are
% included. No dT/dt-bin filter is imposed here. eventT has already restricted
% the reorientation onset to the focused ROI.

    directions=["Coldward","Transverse","Warmward"];

    plateProbT=table();
    summaryT=table();

    % -------------------------------------------------------------
    % First build plate-level conditional transition probabilities.
    % -------------------------------------------------------------
    for gg=1:2
        genotype=string(compareNames(gg));

        G=eventT(string(eventT.Strain)==genotype & ...
            ismember(string(eventT.IncomingDirection),directions) & ...
            ismember(string(eventT.OutgoingDirection),directions),:);

        if isempty(G)
            continue
        end

        plates=unique(string(G.PlateID),'stable');

        for pp=1:numel(plates)
            plate=plates(pp);
            P=G(string(G.PlateID)==plate,:);

            for ii=1:3
                incoming=directions(ii);
                R=P(string(P.IncomingDirection)==incoming,:);

                nRow=height(R);
                if nRow<1
                    continue
                end

                for jj=1:3
                    outgoing=directions(jj);
                    nCell=sum(string(R.OutgoingDirection)==outgoing);
                    prob=nCell/nRow;

                    Genotype=genotype;
                    PlateID=plate;
                    IncomingDirection=incoming;
                    OutgoingDirection=outgoing;
                    TransitionCount=nCell;
                    IncomingRowTotal=nRow;
                    Probability=prob;

                    row=table(Genotype,PlateID,IncomingDirection, ...
                        OutgoingDirection,TransitionCount,IncomingRowTotal, ...
                        Probability);

                    plateProbT=local_append_table(plateProbT,row);
                end
            end
        end
    end

    % -------------------------------------------------------------
    % Summarize equally across independent plates.
    % -------------------------------------------------------------
    for gg=1:2
        genotype=string(compareNames(gg));

        for ii=1:3
            incoming=directions(ii);

            % Same contributing plates are used for all three outgoing cells
            % in one row because plateProbT only creates a row when that plate
            % had >=1 valid transition from the incoming state.
            rowMask=string(plateProbT.Genotype)==genotype & ...
                string(plateProbT.IncomingDirection)==incoming;

            rowData=plateProbT(rowMask,:);
            if isempty(rowData)
                continue
            end

            contributingPlates=unique(string(rowData.PlateID));
            nPlates=numel(contributingPlates);

            for jj=1:3
                outgoing=directions(jj);
                D=rowData(string(rowData.OutgoingDirection)==outgoing,:);

                vals=D.Probability;
                vals=vals(isfinite(vals));

                meanP=mean(vals,'omitnan');
                if numel(vals)>=2
                    semP=std(vals,0,'omitnan')/sqrt(numel(vals));
                else
                    semP=NaN;
                end

                Genotype=genotype;
                ROI_Low_C=roiC(1);
                ROI_High_C=roiC(2);
                IncomingDirection=incoming;
                OutgoingDirection=outgoing;
                MeanProbabilityAcrossPlates=meanP;
                SEMAcrossPlates=semP;
                NPlatesContributing=nPlates;
                TotalTransitionEvents=sum(D.TransitionCount,'omitnan');

                row=table(Genotype,ROI_Low_C,ROI_High_C, ...
                    IncomingDirection,OutgoingDirection, ...
                    MeanProbabilityAcrossPlates,SEMAcrossPlates, ...
                    NPlatesContributing,TotalTransitionEvents);

                summaryT=local_append_table(summaryT,row);
            end
        end
    end

    % -------------------------------------------------------------
    % Plot ONLY the two genotypes relevant to this focused ROI.
    % -------------------------------------------------------------
    fig=figure('Color','w','Position',[90 120 1100 520], ...
        'Name',sprintf('ROI_TransitionMatrices_%.0f_%.0fC',roiC));

    tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    for gg=1:2
        genotype=string(compareNames(gg));
        ax=nexttile(tl,gg);

        M=nan(3,3);
        Nrow=zeros(3,1);

        for ii=1:3
            incoming=directions(ii);

            for jj=1:3
                outgoing=directions(jj);

                q=string(summaryT.Genotype)==genotype & ...
                    string(summaryT.IncomingDirection)==incoming & ...
                    string(summaryT.OutgoingDirection)==outgoing;

                if any(q)
                    M(ii,jj)=summaryT.MeanProbabilityAcrossPlates(find(q,1));
                    Nrow(ii)=summaryT.NPlatesContributing(find(q,1));
                end
            end
        end

        imagesc(ax,M,[0 1]);
        axis(ax,'square');
        colormap(ax,parula(256));
        cb=colorbar(ax);
        cb.Label.String='Probability';
        cb.Ticks=0:0.2:1;

        xticks(ax,1:3);
        yticks(ax,1:3);
        xticklabels(ax,directions);
        yticklabels(ax,directions);
        xtickangle(ax,28);

        xlabel(ax,'Outgoing direction');
        ylabel(ax,'Incoming direction');

        sIdx=find(string(strainNames)==genotype,1);
        if isempty(sIdx)
            titleColor=[0 0 0];
        else
            titleColor=strainColors(sIdx,:);
        end

        title(ax,char(genotype), ...
            'Color',titleColor, ...
            'FontWeight','bold', ...
            'Interpreter','none');

        % Cell probability values.
        for ii=1:3
            for jj=1:3
                if isfinite(M(ii,jj))
                    text(ax,jj,ii,sprintf('%.2f',M(ii,jj)), ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontWeight','bold', ...
                        'FontSize',12, ...
                        'Color',[0.08 0.08 0.08]);
                else
                    text(ax,jj,ii,'NA', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize',10, ...
                        'Color',[0.35 0.35 0.35]);
                end
            end
        end

        % Show how many independent plates contribute to each conditional row.
        rowLabelText=strings(3,1);
        for ii=1:3
            rowLabelText(ii)=sprintf('%s  (n=%d plates)',directions(ii),Nrow(ii));
        end
        yticklabels(ax,rowLabelText);

        box(ax,'on');
        set(ax,'TickDir','out','FontSize',10.5);
    end

    title(tl,sprintf([ ...
        'ROI-specific reorientation transition matrices | %.0f-%.0f °C\n' ...
        'each row = mean across plates of P(outgoing direction | incoming direction, reorientation)'], ...
        roiC(1),roiC(2)), ...
        'FontWeight','bold');
end


function [eventT,exposureT] = local_focus_build_reorientation_data( ...
    coordsByStrain,plateIDdata,plateOrder,strainNames,samplingRate, ...
    pixRange,tempRange,roiC,dTdtWindow_s,headingWindow_s,dTdtEdges, ...
    reversalGap_s,clusterMode,coldHalfWidth,warmThreshold)
% Build clustered reorientation events and matched exposure.
%
% One reorientation event = first call of one reversal-call cluster.
%
% PRECEDING dT/dt:
%   mean signed dT/dt during the 2 s immediately BEFORE candidate/event onset.
%
% FOLLOWING dT/dt:
%   mean signed dT/dt during the 2 s immediately AFTER candidate/event onset.
%   This is a reverse-time association and must NOT be interpreted as a
%   causal predictor of the reorientation.
%
% For rate calculations, every valid processed sample in the ROI is treated
% as a candidate onset opportunity. These candidate samples are binned by
% either PRE or FOLLOWING dT/dt, so both rate plots have matched exposure
% denominators.
%
% At the current processed sampling rate of 1 Hz:
%   probability per opportunity = clustered onsets / eligible 1-s samples.

    CperPx=diff(tempRange)/diff(pixRange);
    winDT=max(1,round(dTdtWindow_s*samplingRate));
    winH=max(1,round(headingWindow_s*samplingRate));
    gapFrames=max(0,round(reversalGap_s*samplingRate));

    nBins=numel(dTdtEdges)-1;
    binCenters=(dTdtEdges(1:end-1)+dTdtEdges(2:end))/2;

    eventT=table();
    exposureT=table();

    for s=1:numel(strainNames)
        coords=coordsByStrain{s};
        trackPlate=string(plateIDdata{s}(:));
        plates=string(plateOrder{s}(:));

        for p=1:numel(plates)
            trIdx=find(trackPlate==plates(p));

            % Total valid tracked animal-time in the focused ROI.
            % This denominator is used for the standalone OVERALL ROI
            % reorientation rate and does NOT require a valid heading or dT/dt.
            expROIAll=0;

            % State-qualified exposure used by the directional analyses below.
            expAll=0;
            expDir=zeros(3,1);
            expDirBin=zeros(3,nBins);
            expColdBin=zeros(nBins,1);
            expPreAllBin=zeros(nBins,1);
            expPostAllBin=zeros(nBins,1);

            for ii=1:numel(trIdx)
                w=trIdx(ii);
                if w>numel(coords) || isempty(coords(w).x_smooth),continue;end

                xcell=coords(w).x_smooth;
                ycell=coords(w).y_smooth;
                rcell=coords(w).rev_index;
                if ~iscell(xcell),xcell={xcell};end
                if ~iscell(ycell),ycell={ycell};end
                if ~iscell(rcell),rcell={rcell};end

                for seg=1:numel(xcell)
                    if seg>numel(ycell),continue;end
                    x=xcell{seg}; y=ycell{seg};
                    if isempty(x)||isempty(y),continue;end

                    x=x(:); y=y(:);
                    N=min(numel(x),numel(y));
                    x=x(1:N); y=y(1:N);

                    temp=tempRange(1)+(x-pixRange(1))*CperPx;

                    % OVERALL ROI TRACKED EXPOSURE:
                    % every finite X/Y sample inside the focused ROI contributes.
                    % At 1-Hz processed sampling, one valid sample ~= 1 animal-s.
                    validROI=isfinite(x) & isfinite(y) & ...
                        isfinite(temp) & temp>=roiC(1) & temp<roiC(2);
                    expROIAll=expROIAll+sum(validROI);

                    % The sensory/heading analyses below require enough samples
                    % for their requested pre/post windows.
                    if N<max(winDT,winH)+3,continue;end

                    % stepDTdt(k) = dT/dt for interval k-1 -> k.
                    stepDTdt=nan(N,1);
                    pairGood=isfinite(x(1:end-1)) & isfinite(x(2:end)) & ...
                        isfinite(y(1:end-1)) & isfinite(y(2:end));
                    step=diff(x)*CperPx*samplingRate;
                    tmp=nan(N-1,1);
                    tmp(pairGood)=step(pairGood);
                    stepDTdt(2:end)=tmp;

                    preceding=nan(N,1);
                    following=nan(N,1);

                    for k=(winDT+1):N
                        vv=stepDTdt(k-winDT+1:k);
                        if numel(vv)==winDT && all(isfinite(vv))
                            preceding(k)=mean(vv);
                        end
                    end

                    for k=1:(N-winDT)
                        vv=stepDTdt(k+1:k+winDT);
                        if numel(vv)==winDT && all(isfinite(vv))
                            following(k)=mean(vv);
                        end
                    end

                    h=local_heading_segment_for_reversal(x,y);

                    % All-direction exposure for PRE and POST/future dT/dt.
                    for k=1:N
                        if ~local_focus_temp_in_roi(temp(k),roiC),continue;end

                        if isfinite(preceding(k))
                            bb=discretize(preceding(k),dTdtEdges);
                            if isfinite(bb),expPreAllBin(bb)=expPreAllBin(bb)+1;end
                        end

                        if isfinite(following(k))
                            bb=discretize(following(k),dTdtEdges);
                            if isfinite(bb),expPostAllBin(bb)=expPostAllBin(bb)+1;end
                        end
                    end

                    % Incoming-direction exposure for coldward interruption.
                    for k=(max(winDT,winH)+1):N
                        if ~local_focus_temp_in_roi(temp(k),roiC),continue;end
                        if ~isfinite(preceding(k)),continue;end

                        h1=k-winH; h2=k-1;
                        if h1<1 || h2>numel(h),continue;end
                        hh=local_circ_mean_deg(h(h1:h2));
                        if ~isfinite(hh),continue;end

                        d=local_focus_direction(hh,coldHalfWidth,warmThreshold);
                        expAll=expAll+1;

                        bb=discretize(preceding(k),dTdtEdges);

                        if d=="Coldward"
                            expDir(1)=expDir(1)+1;
                            if isfinite(bb)
                                expColdBin(bb)=expColdBin(bb)+1;
                                expDirBin(1,bb)=expDirBin(1,bb)+1;
                            end
                        elseif d=="Transverse"
                            expDir(2)=expDir(2)+1;
                            if isfinite(bb),expDirBin(2,bb)=expDirBin(2,bb)+1;end
                        else
                            expDir(3)=expDir(3)+1;
                            if isfinite(bb),expDirBin(3,bb)=expDirBin(3,bb)+1;end
                        end
                    end

                    % Clustered reorientation events.
                    rr=[];
                    if seg<=numel(rcell) && ~isempty(rcell{seg})
                        rr=rcell{seg}(:);
                        rr=sort(unique(round(rr(isfinite(rr)&rr>=1&rr<=N))));
                    end
                    if isempty(rr),continue;end

                    cl=local_cluster_reversal_calls(rr,gapFrames,string(clusterMode));

                    for cc=1:size(cl,1)
                        onset=cl(cc,1);
                        episodeEnd=cl(cc,2);

                        if onset<1 || onset>N || ...
                                ~local_focus_temp_in_roi(temp(onset),roiC)
                            continue
                        end

                        hin=NaN;
                        incoming="Missing";
                        pre1=onset-winH;
                        pre2=onset-1;
                        if pre1>=1 && pre2<=numel(h)
                            hin=local_circ_mean_deg(h(pre1:pre2));
                            if isfinite(hin)
                                incoming=local_focus_direction( ...
                                    hin,coldHalfWidth,warmThreshold);
                            end
                        end

                        hout=NaN;
                        outgoing="Missing";
                        post1=episodeEnd+1;
                        post2=episodeEnd+winH;
                        if post1>=1 && post2<=numel(h)
                            hout=local_circ_mean_deg(h(post1:post2));
                            if isfinite(hout)
                                outgoing=local_focus_direction( ...
                                    hout,coldHalfWidth,warmThreshold);
                            end
                        end

                        preVal=preceding(onset);
                        postVal=following(onset);

                        preBin=NaN; postBin=NaN;
                        if isfinite(preVal),preBin=discretize(preVal,dTdtEdges);end
                        if isfinite(postVal),postBin=discretize(postVal,dTdtEdges);end

                        Strain=string(strainNames{s});
                        PlateID=plates(p);
                        TrackID=w;
                        Segment=seg;
                        OnsetIndex=onset;
                        EndIndex=episodeEnd;
                        OnsetTemperature_C=temp(onset);
                        Preceding_dTdt_C_per_s=preVal;
                        Following_dTdt_C_per_s=postVal;
                        Preceding_dTdtBinIndex=preBin;
                        Following_dTdtBinIndex=postBin;
                        dTdtBinIndex=preBin; % compatibility alias
                        IncomingHeading_deg=hin;
                        IncomingDirection=incoming;
                        OutgoingHeading_deg=hout;
                        OutgoingDirection=outgoing;
                        ClusterSize=cl(cc,3);

                        row=table(Strain,PlateID,TrackID,Segment,OnsetIndex,EndIndex, ...
                            OnsetTemperature_C,Preceding_dTdt_C_per_s, ...
                            Following_dTdt_C_per_s,Preceding_dTdtBinIndex, ...
                            Following_dTdtBinIndex,dTdtBinIndex,IncomingHeading_deg, ...
                            IncomingDirection,OutgoingHeading_deg,OutgoingDirection, ...
                            ClusterSize);

                        eventT=local_append_table(eventT,row);
                    end
                end
            end

            % Per-bin exposure rows.
            dirLabels=["Coldward","Transverse","Warmward"];
            for bb=1:nBins
                Strain=string(strainNames{s});
                PlateID=plates(p);
                dTdtBinIndex=bb;
                dTdtBinCenter_C_per_s=binCenters(bb);

                for dd=1:3
                    Direction=dirLabels(dd);
                    ExposureFrames=expDirBin(dd,bb);
                    Exposure_s=ExposureFrames/samplingRate;
                    Exposure_min=Exposure_s/60;
                    exposureT=local_append_table(exposureT,table( ...
                        Strain,PlateID,Direction,dTdtBinIndex, ...
                        dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));
                end

                Direction="All_PreBin";
                ExposureFrames=expPreAllBin(bb);
                Exposure_s=ExposureFrames/samplingRate;
                Exposure_min=Exposure_s/60;
                exposureT=local_append_table(exposureT,table( ...
                    Strain,PlateID,Direction,dTdtBinIndex, ...
                    dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));

                Direction="All_PostBin";
                ExposureFrames=expPostAllBin(bb);
                Exposure_s=ExposureFrames/samplingRate;
                Exposure_min=Exposure_s/60;
                exposureT=local_append_table(exposureT,table( ...
                    Strain,PlateID,Direction,dTdtBinIndex, ...
                    dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));
            end

            % Aggregate exposure sentinels.
            Strain=string(strainNames{s});
            PlateID=plates(p);

            % TRUE overall ROI exposure: all valid tracked animal-time in ROI,
            % independent of heading or dT/dt availability.
            Direction="ROI_AllTracked";
            dTdtBinIndex=-2;
            dTdtBinCenter_C_per_s=NaN;
            ExposureFrames=expROIAll;
            Exposure_s=ExposureFrames/samplingRate;
            Exposure_min=Exposure_s/60;
            exposureT=local_append_table(exposureT,table( ...
                Strain,PlateID,Direction,dTdtBinIndex, ...
                dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));

            % State-qualified aggregate exposure retained for legacy/directional
            % analyses that require a valid preceding dT/dt and incoming heading.
            Direction="All";
            dTdtBinIndex=0;
            dTdtBinCenter_C_per_s=NaN;
            ExposureFrames=expAll;
            Exposure_s=ExposureFrames/samplingRate;
            Exposure_min=Exposure_s/60;
            exposureT=local_append_table(exposureT,table( ...
                Strain,PlateID,Direction,dTdtBinIndex, ...
                dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));

            dirLabels=["Coldward","Transverse","Warmward"];
            for dd=1:3
                Direction=dirLabels(dd);
                dTdtBinIndex=-1;
                dTdtBinCenter_C_per_s=NaN;
                ExposureFrames=expDir(dd);
                Exposure_s=ExposureFrames/samplingRate;
                Exposure_min=Exposure_s/60;
                exposureT=local_append_table(exposureT,table( ...
                    Strain,PlateID,Direction,dTdtBinIndex, ...
                    dTdtBinCenter_C_per_s,ExposureFrames,Exposure_s,Exposure_min));
            end
        end
    end
end


function T = local_focus_reorientation_rate_by_dTdt( ...
    eventT,exposureT,plateOrder,strainNames,dTdtEdges,minExposure_s,windowMode)
% All-direction clustered reorientation rate versus PRE or FOLLOWING dT/dt.

    windowMode=string(windowMode);

    if windowMode=="Pre"
        directionTag="All_PreBin";
        eventBinVar='Preceding_dTdtBinIndex';
    elseif windowMode=="Post"
        directionTag="All_PostBin";
        eventBinVar='Following_dTdtBinIndex';
    else
        error('windowMode must be "Pre" or "Post".');
    end

    nBins=numel(dTdtEdges)-1;
    centers=(dTdtEdges(1:end-1)+dTdtEdges(2:end))/2;
    T=table();

    for s=1:numel(strainNames)
        plates=string(plateOrder{s}(:));

        for p=1:numel(plates)
            plate=plates(p);

            e=eventT(string(eventT.Strain)==string(strainNames{s}) & ...
                string(eventT.PlateID)==plate,:);

            for bb=1:nBins
                x=exposureT( ...
                    string(exposureT.Strain)==string(strainNames{s}) & ...
                    string(exposureT.PlateID)==plate & ...
                    string(exposureT.Direction)==directionTag & ...
                    exposureT.dTdtBinIndex==bb,:);

                if isempty(x),continue;end

                eb=e(isfinite(e.(eventBinVar)) & e.(eventBinVar)==bb,:);

                Strain=string(strainNames{s});
                PlateID=plate;
                Window=windowMode;
                BinIndex=bb;
                BinCenter=centers(bb);
                Count=height(eb);
                ExposureFrames=x.ExposureFrames(1);
                Exposure_s=x.Exposure_s(1);
                Exposure_min=x.Exposure_min(1);

                Included=Exposure_s>=minExposure_s && ...
                    ExposureFrames>0 && Exposure_min>0;

                Rate=NaN;
                ProbabilityPerOpportunity=NaN;
                ProbabilityPercentPerOpportunity=NaN;

                if Included
                    Rate=Count/Exposure_min;
                    ProbabilityPerOpportunity=Count/ExposureFrames;
                    ProbabilityPercentPerOpportunity=100*ProbabilityPerOpportunity;
                end

                row=table(Strain,PlateID,Window,BinIndex,BinCenter,Count, ...
                    ExposureFrames,Exposure_s,Exposure_min,Included,Rate, ...
                    ProbabilityPerOpportunity,ProbabilityPercentPerOpportunity);

                T=local_append_table(T,row);
            end
        end
    end
end


function T = local_focus_true_direction_exit_rate_table( ...
    eventT,exposureT,plateOrder,strainNames,minExposure_s,incomingDirection)
% STRICT DIRECTION-EXIT RATE
%
% Numerator:
%   clustered reorientation onsets that begin in incomingDirection and end
%   OUTSIDE that direction.
%
% Denominator:
%   tracked animal-minutes spent in incomingDirection.
%
% For Warmward:
%   true warmward-exit rate =
%       #(Warmward -> Transverse OR Coldward)
%       -------------------------------------
%       warmward animal-minutes
%
% For Coldward:
%   true coldward-exit rate =
%       #(Coldward -> Transverse OR Warmward)
%       -------------------------------------
%       coldward animal-minutes

    incomingDirection=string(incomingDirection);

    if incomingDirection=="Warmward"
        allowedOutgoing=["Coldward","Transverse"];
    elseif incomingDirection=="Coldward"
        allowedOutgoing=["Transverse","Warmward"];
    elseif incomingDirection=="Transverse"
        allowedOutgoing=["Coldward","Warmward"];
    else
        error('incomingDirection must be Coldward, Transverse, or Warmward.');
    end

    T=table();

    for s=1:numel(strainNames)
        plates=string(plateOrder{s}(:));

        for p=1:numel(plates)
            plate=plates(p);

            e=eventT( ...
                string(eventT.Strain)==string(strainNames{s}) & ...
                string(eventT.PlateID)==plate & ...
                string(eventT.IncomingDirection)==incomingDirection & ...
                ismember(string(eventT.OutgoingDirection),allowedOutgoing),:);

            x=exposureT( ...
                string(exposureT.Strain)==string(strainNames{s}) & ...
                string(exposureT.PlateID)==plate & ...
                string(exposureT.Direction)==incomingDirection & ...
                exposureT.dTdtBinIndex==-1,:);

            if isempty(x),continue;end

            Strain=string(strainNames{s});
            PlateID=plate;
            Direction=incomingDirection;
            ReversalCount=height(e);
            TrackedTime_s=x.Exposure_s(1);
            TrackedTime_min=x.Exposure_min(1);

            IncludedForRate=TrackedTime_s>=minExposure_s && ...
                isfinite(TrackedTime_min) && TrackedTime_min>0;

            ReversalRate_per_tracked_min=NaN;
            if IncludedForRate
                ReversalRate_per_tracked_min=ReversalCount/TrackedTime_min;
            end

            row=table(Strain,PlateID,Direction,ReversalCount,TrackedTime_s, ...
                TrackedTime_min,ReversalRate_per_tracked_min,IncludedForRate);

            T=local_append_table(T,row);
        end
    end
end


function T = local_focus_build_clustered_complete_runs( ...
    coordsByStrain,plateIDdata,strainNames,samplingRate,pixRange,tempRange,roiC, ...
    reversalGap_s,clusterMode,coldHalfWidth,warmThreshold,minRun_s,runSpeedCutoff)
% Complete inter-reorientation runs using CLUSTERED reorientation episodes.
%
% A complete run begins after the LAST call of one reorientation cluster and
% ends before the FIRST call of the next cluster. This avoids counting the
% multi-call reorientation episode itself as part of the run.
%
% Only runs STARTING inside roiC are returned.
%
% CrossedColderROIBoundary:
%   true if the run reaches a temperature below roiC(1) before the next
%   reorientation. This provides a boundary-defined "successful progress"
%   measure for eat-4 and a "continued through the goal ROI" measure for unc-17.

    CperPx=diff(tempRange)/diff(pixRange);
    gapFrames=max(0,round(reversalGap_s*samplingRate));
    T=table();

    for s=1:numel(strainNames)
        coords=coordsByStrain{s};
        plates=string(plateIDdata{s}(:));

        for w=1:min(numel(coords),numel(plates))
            if isempty(coords(w).x_smooth),continue;end
            xcell=coords(w).x_smooth; ycell=coords(w).y_smooth; rcell=coords(w).rev_index;
            if ~iscell(xcell),xcell={xcell};end
            if ~iscell(ycell),ycell={ycell};end
            if ~iscell(rcell),rcell={rcell};end

            for seg=1:min(numel(xcell),numel(ycell))
                x=xcell{seg}; y=ycell{seg};
                if isempty(x)||isempty(y),continue;end
                x=x(:);y=y(:);N=min(numel(x),numel(y));x=x(1:N);y=y(1:N);

                rr=[];
                if seg<=numel(rcell)&&~isempty(rcell{seg})
                    rr=rcell{seg}(:);
                    rr=sort(unique(round(rr(isfinite(rr)&rr>=1&rr<=N))));
                end
                if numel(rr)<2,continue;end
                cl=local_cluster_reversal_calls(rr,gapFrames,string(clusterMode));
                if size(cl,1)<2,continue;end

                for c=1:(size(cl,1)-1)
                    rs=cl(c,2)+1;
                    re=cl(c+1,1)-1;
                    if rs<1||re>N||re<=rs,continue;end

                    xr=x(rs:re); yr=y(rs:re);
                    if any(~isfinite(xr)|~isfinite(yr)),continue;end

                    duration=(re-rs)/samplingRate;
                    if duration<minRun_s,continue;end

                    dxy=hypot(diff(xr),diff(yr));
                    dist_mm=sum(dxy,'omitnan')/(diff(pixRange)/100);
                    meanSpeed=dist_mm/duration;
                    if ~isfinite(meanSpeed)||meanSpeed<runSpeedCutoff,continue;end

                    startT=tempRange(1)+(xr(1)-pixRange(1))*CperPx;
                    endT=tempRange(1)+(xr(end)-pixRange(1))*CperPx;
                    if ~local_focus_temp_in_roi(startT,roiC),continue;end

                    [hh,~]=calculateRelativeHeadingandDistance(xr,yr,[-1 0]);
                    meanHeading=local_circ_mean_deg(hh);
                    orientation=local_focus_direction(meanHeading,coldHalfWidth,warmThreshold);

                    tt=tempRange(1)+(xr-pixRange(1))*CperPx;
                    crossed=any(tt<roiC(1));

                    Strain=string(strainNames{s});
                    PlateID=plates(w);
                    TrackID=w;
                    Segment=seg;
                    RunIndex=c;
                    RunStart=rs;
                    RunEnd=re;
                    StartTemperature_C=startT;
                    EndTemperature_C=endT;
                    SignedThermalDisplacement_C=endT-startT;
                    RunDuration_s=duration;
                    RunDistance_mm=dist_mm;
                    RunMeanSpeed_mm_s=meanSpeed;
                    RunHeading_deg=meanHeading;
                    Orientation=orientation;
                    CrossedColderROIBoundary=double(crossed);

                    row=table(Strain,PlateID,TrackID,Segment,RunIndex,RunStart,RunEnd, ...
                        StartTemperature_C,EndTemperature_C,SignedThermalDisplacement_C, ...
                        RunDuration_s,RunDistance_mm,RunMeanSpeed_mm_s,RunHeading_deg, ...
                        Orientation,CrossedColderROIBoundary);
                    T=local_append_table(T,row);
                end
            end
        end
    end
end


function [fig,selectedT] = local_focus_plot_top_directional_runs( ...
    runT,coordsByStrain,compareNames,strainNames,samplingRate, ...
    pixRange,tempRange,roiC,topN,coldColor,warmColor)
% Visualize the TOP-X LONGEST complete Coldward and Warmward runs.
%
% Selection:
%   Within each genotype x direction, sort eligible complete runs by
%   RunDuration_s from longest to shortest and plot the first topN.
%
% Visualization:
%   Every selected run is translated so its first valid coordinate is (0,0).
%   The trajectories are NOT rotated, reflected, stretched, or normalized.
%   Therefore the plotted path preserves the animal's REAL X/Y displacement.
%
% Coordinate convention:
%   X increases toward warmer temperatures.
%   Therefore coldward runs generally extend toward NEGATIVE relative X,
%   whereas warmward runs generally extend toward POSITIVE relative X.
%
% Units:
%   Relative X and Y displacement are shown in mm.
%
% Layout:
%   row 1 = genotype 1
%   row 2 = genotype 2
%   col 1 = Coldward runs
%   col 2 = Warmward runs
%
% NOTE:
%   runT already contains only complete qualifying runs that START inside
%   the focused ROI and satisfy the focused duration/speed requirements.

    %#ok<INUSD> samplingRate tempRange roiC

    fig=figure('Color','w','Position',[70 60 1320 800], ...
        'Name',sprintf('TopDirectionalRuns_%.0f_%.0fC',roiC));

    directions=["Coldward","Warmward"];
    dirColors=[coldColor;warmColor];

    % Physical calibration used throughout this assay:
    % 900 px across the 15-25 C gradient = 100 mm arena width.
    pixelsPerMm=diff(pixRange)/100;

    % Store the selected metadata and actual centered trajectories before
    % plotting so ALL four panels can share identical X/Y limits.
    selectedT=table();
    traceX=cell(2,2);
    traceY=cell(2,2);
    traceDuration=cell(2,2);

    globalMaxAbsX=0;
    globalMaxAbsY=0;

    for gg=1:2
        genotype=string(compareNames(gg));
        sIdx=find(string(strainNames)==genotype,1);

        for dd=1:2
            direction=directions(dd);
            traceX{gg,dd}={};
            traceY{gg,dd}={};
            traceDuration{gg,dd}=[];

            if isempty(sIdx) || isempty(runT)
                continue
            end

            D=runT(string(runT.Strain)==genotype & ...
                string(runT.Orientation)==direction,:);

            if isempty(D)
                continue
            end

            % TOP = LONGEST DURATION.
            D=sortrows(D,'RunDuration_s','descend');

            nShow=min(max(round(topN),1),height(D));
            D=D(1:nShow,:);
            D.PlotRank=(1:nShow)';
            D.SelectionDefinition=repmat( ...
                "ranked by longest RunDuration_s",nShow,1);

            if isempty(selectedT)
                selectedT=D;
            else
                selectedT=[selectedT;D]; %#ok<AGROW>
            end

            coords=coordsByStrain{sIdx};

            for rr=1:nShow
                w=round(D.TrackID(rr));
                seg=round(D.Segment(rr));
                rs=round(D.RunStart(rr));
                re=round(D.RunEnd(rr));

                if w<1 || w>numel(coords)
                    continue
                end
                if ~isfield(coords,'x_smooth') || ~isfield(coords,'y_smooth') || ...
                        isempty(coords(w).x_smooth) || isempty(coords(w).y_smooth)
                    continue
                end

                xcell=coords(w).x_smooth;
                ycell=coords(w).y_smooth;
                if ~iscell(xcell),xcell={xcell};end
                if ~iscell(ycell),ycell={ycell};end

                if seg<1 || seg>numel(xcell) || seg>numel(ycell) || ...
                        isempty(xcell{seg}) || isempty(ycell{seg})
                    continue
                end

                x=xcell{seg}(:);
                y=ycell{seg}(:);
                N=min(numel(x),numel(y));
                x=x(1:N);
                y=y(1:N);

                if rs<1 || re>N || re<rs
                    continue
                end

                xr=x(rs:re);
                yr=y(rs:re);

                good=isfinite(xr) & isfinite(yr);
                if sum(good)<2
                    continue
                end

                % Preserve the actual trajectory geometry, but translate the
                % first valid point to a shared origin.
                firstGood=find(good,1,'first');

                dx_mm=(xr-xr(firstGood))/pixelsPerMm;
                dy_mm=(yr-yr(firstGood))/pixelsPerMm;

                dx_mm=dx_mm(good);
                dy_mm=dy_mm(good);

                traceX{gg,dd}{end+1,1}=dx_mm; %#ok<AGROW>
                traceY{gg,dd}{end+1,1}=dy_mm; %#ok<AGROW>
                traceDuration{gg,dd}(end+1,1)=D.RunDuration_s(rr); %#ok<AGROW>

                globalMaxAbsX=max(globalMaxAbsX,max(abs(dx_mm),[],'omitnan'));
                globalMaxAbsY=max(globalMaxAbsY,max(abs(dy_mm),[],'omitnan'));
            end
        end
    end

    % Common axes across ALL four panels so displacement magnitudes are
    % directly comparable between genotypes and directions.
    if ~isfinite(globalMaxAbsX) || globalMaxAbsX<=0
        globalMaxAbsX=1;
    end
    if ~isfinite(globalMaxAbsY) || globalMaxAbsY<=0
        globalMaxAbsY=1;
    end

    xLim=1.08*globalMaxAbsX*[-1 1];
    yLim=1.08*globalMaxAbsY*[-1 1];

    tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    for gg=1:2
        genotype=string(compareNames(gg));

        for dd=1:2
            direction=directions(dd);
            ax=nexttile(tl,(gg-1)*2+dd);
            hold(ax,'on');

            Xs=traceX{gg,dd};
            Ys=traceY{gg,dd};

            % Common point of origin.
            xline(ax,0,':','Color',[0.55 0.55 0.55], ...
                'LineWidth',0.9,'HandleVisibility','off');
            yline(ax,0,':','Color',[0.55 0.55 0.55], ...
                'LineWidth',0.9,'HandleVisibility','off');

            nPlotted=numel(Xs);

            for rr=1:nPlotted
                plot(ax,Xs{rr},Ys{rr},'-', ...
                    'Color',dirColors(dd,:), ...
                    'LineWidth',1.45, ...
                    'HandleVisibility','off');
            end

            % Mark the common starting point on top of all traces.
            scatter(ax,0,0,28,'o', ...
                'MarkerFaceColor','w', ...
                'MarkerEdgeColor',[0.20 0.20 0.20], ...
                'LineWidth',0.9, ...
                'HandleVisibility','off');

            xlim(ax,xLim);
            ylim(ax,yLim);
            axis(ax,'equal');

            xlabel(ax,'Relative X displacement (mm)');
            ylabel(ax,'Relative Y displacement (mm)');

            if nPlotted>0
                medDur=median(traceDuration{gg,dd},'omitnan');
                maxDur=max(traceDuration{gg,dd},[],'omitnan');
                title(ax,sprintf('%s | %s | longest %d\nmedian %.1f s | max %.1f s', ...
                    genotype,direction,nPlotted,medDur,maxDur), ...
                    'Interpreter','none','FontWeight','normal');
            else
                title(ax,sprintf('%s | %s | no eligible runs', ...
                    genotype,direction), ...
                    'Interpreter','none','FontWeight','normal');
            end

            box(ax,'off');
            set(ax,'TickDir','out');
            ax.XGrid='on';
            ax.YGrid='on';
            ax.GridAlpha=0.12;
        end
    end

    title(tl,sprintf([ ...
        'Focused ROI %.0f-%.0f °C | longest complete directional runs\n' ...
        'all trajectories translated to the same origin; no rotation or rescaling'], ...
        roiC(1),roiC(2)), ...
        'FontWeight','normal');
end


function R = local_focus_fit_rate_glmm(T,compareNames)
% Two-group Poisson rate model with exposure offset.
%
% Each ROI contributes ONE count/exposure observation per independent plate,
% so no repeated plate term is needed here. Plate is already the observation
% unit. A Poisson GLM is therefore preferable to adding a random intercept
% with only one observation per plate.

    if isempty(T) || ~all(ismember({'Strain','IncludedForRate','TrackedTime_min', ...
            'ReversalCount'},T.Properties.VariableNames))
        R=local_focus_empty_fit();
        return
    end

    keep=ismember(string(T.Strain),string(compareNames)) & ...
        T.IncludedForRate & T.TrackedTime_min>0 & isfinite(T.ReversalCount);
    D=T(keep,:);

    if isempty(D) || numel(unique(string(D.Strain)))<2
        R=local_focus_empty_fit();
        return
    end

    D.Genotype=categorical(string(D.Strain),string(compareNames),string(compareNames));
    D.LogExposure=log(D.TrackedTime_min);

    M=fitglm(D,'ReversalCount ~ Genotype', ...
        'Distribution','poisson','Link','log','Offset',D.LogExposure);

    [b,se,p]=local_focus_genotype_coefficient(M);

    R=struct();
    R.Model=M;
    R.ModelLabel="Poisson GLM; one independent observation per plate";
    R.Estimate=exp(b);
    R.CI95Low=exp(b-1.959963984540054*se);
    R.CI95High=exp(b+1.959963984540054*se);
    R.P=p;
    R.EffectLabel="Rate ratio mutant / N2";
    R.AnalysisData=D;
end


function R = local_focus_fit_rate_tuning_glmm( ...
    T,compareNames,dTdtEdges,minPlatesPerGenotypeBin)
% Poisson GLMM for binned reorientation onset counts with exposure offset.
% Used for both all-direction sensory-rate tables and coldward-interruption
% tables.
%
% A dT/dt bin enters inference only if BOTH genotypes contribute at least
% minPlatesPerGenotypeBin independent plates after exposure filtering.

    if nargin<4 || isempty(minPlatesPerGenotypeBin)
        minPlatesPerGenotypeBin=4;
    end

    if isempty(T) || ~istable(T) || ...
            ~all(ismember({'Strain','Included','Exposure_min','BinIndex','Count'}, ...
            T.Properties.VariableNames))
        R=local_focus_empty_fit();
        return
    end

    D=T(ismember(string(T.Strain),string(compareNames)) & T.Included & ...
        T.Exposure_min>0 & isfinite(T.BinIndex),:);
    if isempty(D),R=local_focus_empty_fit();return;end

    bins=sort(unique(D.BinIndex));
    validBins=[];
    for b=bins(:)'
        q=D.BinIndex==b;
        if numel(unique(string(D.PlateID(q & ...
                string(D.Strain)==string(compareNames(1))))))>=minPlatesPerGenotypeBin && ...
           numel(unique(string(D.PlateID(q & ...
                string(D.Strain)==string(compareNames(2))))))>=minPlatesPerGenotypeBin
            validBins(end+1)=b; %#ok<AGROW>
        end
    end
    D=D(ismember(D.BinIndex,validBins),:);
    if numel(validBins)<2,R=local_focus_empty_fit();return;end

    D.Genotype=categorical(string(D.Strain),string(compareNames),string(compareNames));
    D.BinFactor=categorical(D.BinIndex,validBins,string(validBins));
    D.PlateUID=categorical(string(D.Strain)+"__"+string(D.PlateID));
    D.ObservationID=categorical("obs_"+string((1:height(D))'));
    D.LogExposure=log(D.Exposure_min);

    base='Count ~ Genotype*BinFactor + (1|PlateUID)';
    M0=fitglme(D,base,'Distribution','Poisson','Link','log', ...
        'Offset',D.LogExposure,'FitMethod','Laplace','DummyVarCoding','reference');

    rr=residuals(M0,'ResidualType','Pearson','Conditional',true);
    df=max(sum(isfinite(rr))-numel(fixedEffects(M0)),1);
    dispersion=sum(rr(isfinite(rr)).^2)/df;

    if isfinite(dispersion)&&dispersion>1.5
        M=fitglme(D,[base ' + (1|ObservationID)'], ...
            'Distribution','Poisson','Link','log','Offset',D.LogExposure, ...
            'FitMethod','Laplace','DummyVarCoding','reference');
    else
        M=M0;
    end

    pInt=local_focus_joint_wald(M,["Genotype","BinFactor"]);

    beta=fixedEffects(M);
    V=M.CoefficientCovariance;
    X=designMatrix(M,'Fixed');

    BT=table();
    for b=validBins(:)'
        qc=D.BinIndex==b & string(D.Strain)==string(compareNames(1));
        qm=D.BinIndex==b & string(D.Strain)==string(compareNames(2));
        if ~any(qc)||~any(qm),continue;end
        xc=mean(X(qc,:),1);
        xm=mean(X(qm,:),1);
        c=xm-xc;
        est=c*beta;
        se=sqrt(max(c*V*c',0));
        p=2*(1-normcdf(abs(est/max(se,eps))));
        BinIndex=b;
        BinCenter=mean(dTdtEdges(b:b+1));
        RateRatio=exp(est);
        CI95Low=exp(est-1.959963984540054*se);
        CI95High=exp(est+1.959963984540054*se);
        P=p;
        BT=local_append_table(BT,table(BinIndex,BinCenter,RateRatio,CI95Low,CI95High,P));
    end
    if ~isempty(BT),BT.P_Holm=local_holm_adjust(BT.P);end

    R=struct('Model',M,'InteractionP',pInt,'BinContrasts',BT, ...
        'AnalysisData',D,'Estimate',NaN,'CI95Low',NaN,'CI95High',NaN, ...
        'P',pInt,'EffectLabel',"Genotype x dT/dt-bin interaction");
end


function [b,se,p] = local_focus_genotype_coefficient(M)
% Extract the non-interaction genotype coefficient robustly across MATLAB
% model classes/releases.
%
% fitglm commonly stores coefficient names in Properties.RowNames, whereas
% fitlme/fitglme may expose a table variable called Name.

    C=local_force_table(M.Coefficients);
    nm=local_focus_coefficient_names(C);

    k=find(contains(nm,"Genotype",'IgnoreCase',true) & ...
        ~contains(nm,":"),1);

    if isempty(k)
        warning('Could not identify a Genotype coefficient. Available coefficients: %s', ...
            strjoin(cellstr(nm),', '));
        b=NaN;se=NaN;p=NaN;
        return
    end

    b=C.Estimate(k);
    se=C.SE(k);

    if ismember('pValue',C.Properties.VariableNames)
        p=C.pValue(k);
    elseif ismember('pvalue',C.Properties.VariableNames)
        p=C.pvalue(k);
    elseif ismember('PValue',C.Properties.VariableNames)
        p=C.PValue(k);
    else
        p=2*(1-normcdf(abs(b/max(se,eps))));
    end
end


function p = local_focus_joint_wald(M,tokens)
% Joint Wald test for all fixed coefficients whose names contain ALL tokens.
% Compatible with fitglm, fitlme and fitglme coefficient-table layouts.

    C=local_force_table(M.Coefficients);
    nm=local_focus_coefficient_names(C);

    idx=true(height(C),1);
    for t=1:numel(tokens)
        idx=idx & contains(nm,string(tokens(t)),'IgnoreCase',true);
    end
    idx=find(idx);

    if isempty(idx)
        warning('Joint Wald test could not find requested coefficient terms: %s', ...
            strjoin(cellstr(string(tokens)),', '));
        p=NaN;
        return
    end

    b=C.Estimate(idx);

    try
        V=M.CoefficientCovariance(idx,idx);
    catch
        try
            V=vcov(M);
            V=V(idx,idx);
        catch
            p=NaN;
            warning('Could not obtain coefficient covariance matrix for Wald test.');
            return
        end
    end

    W=b'*(pinv(V))*b;
    df=numel(b);
    p=1-chi2cdf(W,df);
end


function nm = local_focus_coefficient_names(C)
% Return coefficient names from a MATLAB coefficient table.
%
% MATLAB model classes are inconsistent:
%   fitglm          -> names are commonly table RowNames
%   fitlme/fitglme  -> may include a Name variable

    if ~istable(C)
        error('Coefficient output could not be converted to a MATLAB table.');
    end

    vars=string(C.Properties.VariableNames);

    if any(vars=="Name")
        nm=string(C.Name);
    elseif ~isempty(C.Properties.RowNames)
        nm=string(C.Properties.RowNames);
    elseif any(vars=="Term")
        nm=string(C.Term);
    elseif any(vars=="Coefficient")
        nm=string(C.Coefficient);
    else
        nm="Coefficient_"+string((1:height(C))');
    end

    nm=nm(:);
end


function R = local_focus_empty_fit()
    R=struct('Model',[],'FullModel',[],'Estimate',NaN,'CI95Low',NaN, ...
        'CI95High',NaN,'P',NaN,'InteractionP',NaN, ...
        'EffectLabel',"not estimable",'AnalysisData',table(), ...
        'BinContrasts',table(),'CommonSupport',[NaN NaN]);
end


function P = local_focus_subset_plate_rate(T,compareNames)
    P=table();
    if isempty(T) || ~istable(T) || ...
            ~all(ismember({'Strain','PlateID','IncludedForRate', ...
            'ReversalRate_per_tracked_min','ReversalCount'},T.Properties.VariableNames))
        return
    end
    for g=1:numel(compareNames)
        G=T(string(T.Strain)==string(compareNames(g)) & T.IncludedForRate,:);
        for i=1:height(G)
            Strain=string(compareNames(g));PlateID=string(G.PlateID(i));
            Value=G.ReversalRate_per_tracked_min(i);
            N=G.ReversalCount(i);
            P=local_append_table(P,table(Strain,PlateID,Value,N));
        end
    end
end


function local_focus_plate_dot_panel(ax,P,compareNames,strainNames,strainColors, ...
    ttl,ylab,fit,referenceLine)

    hold(ax,'on');
    if isempty(P) || ~istable(P) || ~all(ismember({'Strain','Value'},P.Properties.VariableNames))
        text(ax,0.5,0.5,'Not estimable','Units','normalized', ...
            'HorizontalAlignment','center');
        title(ax,ttl,'FontWeight','normal');
        ylabel(ax,ylab);
        box(ax,'off');
        return
    end
    for g=1:numel(compareNames)
        s=find(string(strainNames)==string(compareNames(g)),1);
        v=P.Value(string(P.Strain)==string(compareNames(g)));
        v=v(isfinite(v));
        if isempty(v),continue;end
        jitter=linspace(-0.06,0.06,numel(v))';
        scatter(ax,g+jitter,v,30,strainColors(s,:),'filled','MarkerFaceAlpha',0.55);
        mu=mean(v,'omitnan');
        se=std(v,0,'omitnan')/sqrt(numel(v));
        errorbar(ax,g,mu,se,'ko','MarkerFaceColor',strainColors(s,:), ...
            'LineWidth',1.4,'CapSize',6);
    end
    xlim(ax,[0.6 2.4]);xticks(ax,1:2);xticklabels(ax,compareNames);
    ylabel(ax,ylab);title(ax,ttl,'FontWeight','normal');
    if isfinite(referenceLine)
        yline(ax,referenceLine,'k--','HandleVisibility','off');
    end
    txt=sprintf('%s = %.3g [%.3g, %.3g]\\np = %s', ...
        fit.EffectLabel,fit.Estimate,fit.CI95Low,fit.CI95High, ...
        local_focus_ptext(fit.P));
    text(ax,0.03,0.97,txt,'Units','normalized','VerticalAlignment','top', ...
        'FontSize',8,'Interpreter','none');
    box(ax,'off');grid(ax,'off');set(ax,'TickDir','out');
end


function tf = local_focus_temp_in_roi(t,roiC)
    tf=isfinite(t) && t>=roiC(1) && t<roiC(2);
end


function d = local_focus_direction(h,coldHalfWidth,warmThreshold)
    a=abs(local_wrap180(h));
    if a<coldHalfWidth
        d="Coldward";
    elseif a<=warmThreshold
        d="Transverse";
    else
        d="Warmward";
    end
end


function a = local_wrap180(a)
    a=mod(a+180,360)-180;
end


function s = local_focus_ptext(p)
    if ~isfinite(p)
        s='n/a';
    elseif p<1e-4
        s=sprintf('%.2g',p);
    else
        s=sprintf('%.4f',p);
    end
end


function local_save_clean_figure(fig,path)
    if isempty(fig)||~ishandle(fig),return;end
    try
        exportgraphics(fig,path,'ContentType','vector');
    catch
        try
            print(fig,char(path),'-dpdf','-painters','-bestfit');
        catch ME
            warning('Could not save figure %s: %s',path,ME.message);
        end
    end
end
