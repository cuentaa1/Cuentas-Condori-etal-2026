%% PREPROCESS XY TRAJECTORIES AND ALIGN THERMAL PLATES
% =========================================================================
% Repository-facing preprocessing script.
%
% Purpose
% -------
% 1. Read tracked worm trajectories from independent assay-plate folders.
% 2. Retain tracks meeting a minimum length requirement.
% 3. Place each track on a common experiment timeline.
% 4. Align each recording to a common X-coordinate reference frame using
%    the midpoint between the two fixed thermostats as a physical landmark.
%    Because the metal slab can shift between experiments, this midpoint can
%    appear at different image X coordinates. Each recording is translated
%    so that the measured thermostat midpoint is placed at X = 1350 pixels.
% 5. Preserve the independent plate identity and original within-plate track
%    identity of every output trajectory.
% 6. Save X, Y, plate ID, track ID, plate name, and alignment metadata.
%
% THERMAL-PLATE ALIGNMENT
% -----------------------
% The midpoint between the two fixed thermostats is measured in the image for
% each independent assay plate. The desired common reference is X = 1350 px.
%
%       x_shift_px = 1350 - measured_thermostat_center_px
%       X_aligned  = X_raw + x_shift_px
%
% Example:
%   measured center = 1315 px -> shift = +35 px
%   measured center = 1380 px -> shift = -30 px
%
% Y coordinates are not shifted.
%
% Each row of the manifest corresponds to ONE independent biological plate.
% =========================================================================

clc
clear
close all

%% ========================================================================
% USER SETTINGS
% ========================================================================

% Root directory containing the raw assay folders listed in the manifest.
% Use a local path on your computer. Do not commit machine-specific paths.
raw_data_root = "/path/to/raw_tracking_data";

% CSV with columns:
%   plate_folder           relative path from raw_data_root to one assay plate
%   thermostat_center_px   measured midpoint between the two fixed thermostats
manifest_file = "plate_alignment_manifest.csv";

% Common X coordinate to which the thermostat midpoint is aligned.
reference_thermostat_center_px = 1350;

% Output.
output_dir = "./processed_data";
dataset_label = "N2";  % examples: "N2", "eat4", "unc17"

% Assay/tracking parameters.
fps = 2;
total_minutes = 15;

% Minimum original track length required for retention.
% At 2 fps, 300 frames = 150 s = 2.5 min.
min_track_frames = 300;

% Raw tracking file expected inside each plate folder.
mat_file_pattern = "*postclean_experiment_1*.mat";

% Optional quality-control plot of aligned X positions.
make_alignment_qc_plot = true;

%% ========================================================================
% READ AND VALIDATE MANIFEST
% ========================================================================

manifest = readtable(manifest_file,'TextType','string');

requiredVariables = ["plate_folder","thermostat_center_px"];
missingVariables = setdiff(requiredVariables,string(manifest.Properties.VariableNames));

if ~isempty(missingVariables)
    error("Manifest is missing required column(s): %s", ...
        strjoin(missingVariables,", "));
end

if isempty(manifest)
    error("Alignment manifest contains no plates.");
end

if any(ismissing(manifest.plate_folder))
    error("Every manifest row must contain a plate_folder.");
end

if any(~isfinite(manifest.thermostat_center_px))
    error("Every manifest row must contain a finite thermostat_center_px value.");
end

% Compute the plate-specific translation directly from the physical landmark.
manifest.x_shift_px = ...
    reference_thermostat_center_px - manifest.thermostat_center_px;

% Resolve plate folders relative to raw_data_root.
plate_folders = strings(height(manifest),1);

for p = 1:height(manifest)
    plate_folders(p) = fullfile(raw_data_root,manifest.plate_folder(p));
end

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% ========================================================================
% EXTRACT, ALIGN, AND POOL TRAJECTORIES
% ========================================================================

[X_data,Y_data,plate_id,track_id,plate_name, ...
    track_metadata,plate_summary] = ...
    collect_xy_from_folders_correctX_plateID_GITHUB( ...
        plate_folders, ...
        manifest.x_shift_px, ...
        min_track_frames, ...
        total_minutes, ...
        fps, ...
        mat_file_pattern);

%% ========================================================================
% SAVE OUTPUTS IN THE FORMAT USED BY DOWNSTREAM FIGURE ANALYSES
% ========================================================================
% Separate MAT files are written because downstream manuscript scripts load:
%   X_<dataset>.mat
%   Y_<dataset>.mat
%   plate_<dataset>.mat
%
% Additional metadata files retain traceability back to the original plate
% and within-plate track.

xVariable = "X_" + dataset_label;
yVariable = "Y_" + dataset_label;
plateVariable = "plate_" + dataset_label;

SX = struct();
SY = struct();
SP = struct();

SX.(xVariable) = X_data;
SY.(yVariable) = Y_data;
SP.(plateVariable) = plate_id;

save(fullfile(output_dir,xVariable + ".mat"),'-struct','SX');
save(fullfile(output_dir,yVariable + ".mat"),'-struct','SY');
save(fullfile(output_dir,plateVariable + ".mat"),'-struct','SP');

writetable(track_metadata, ...
    fullfile(output_dir,dataset_label + "_track_metadata.csv"));

writetable(plate_summary, ...
    fullfile(output_dir,dataset_label + "_plate_alignment_summary.csv"));

% Preserve the exact input manifest used for this processing run.
manifest_used = manifest;
manifest_used.resolved_plate_folder = plate_folders;

writetable(manifest_used, ...
    fullfile(output_dir,dataset_label + "_alignment_manifest_used.csv"));

% Compact combined MATLAB file with all preprocessing outputs.
save(fullfile(output_dir,dataset_label + "_processed_XY_with_metadata.mat"), ...
    'X_data','Y_data','plate_id','track_id','plate_name', ...
    'track_metadata','plate_summary','manifest_used', ...
    'fps','total_minutes','min_track_frames','mat_file_pattern','-v7.3');

%% ========================================================================
% OPTIONAL ALIGNMENT QC
% ========================================================================
% This plot is descriptive only. It shows aligned X-coordinate distributions
% for each independent plate after the supplied X shifts have been applied.

if make_alignment_qc_plot && ~isempty(X_data)

    fig = figure('Color','w', ...
        'Name',dataset_label + " aligned X-position QC", ...
        'Position',[100 100 850 450]);

    ax = axes(fig);
    hold(ax,'on');

    uniquePlates = unique(plate_id(:))';

    for p = uniquePlates
        idx = plate_id == p;

        % Plot a deterministic subset of time points to keep the QC figure
        % responsive without changing the underlying saved data.
        Xp = X_data(:,idx);
        step = max(1,round(size(Xp,1)/500));
        xvals = Xp(1:step:end,:);
        xvals = xvals(isfinite(xvals));

        if isempty(xvals)
            continue
        end

        % Slight vertical jitter is used only for visualization.
        y = p + 0.10.*sin((1:numel(xvals))');
        scatter(ax,xvals,y,4,'filled', ...
            'MarkerFaceAlpha',0.08, ...
            'MarkerEdgeAlpha',0.08);
    end

    xline(ax,reference_thermostat_center_px,'--', ...
        'Thermostat midpoint reference (1350 px)', ...
        'Color',[0.25 0.25 0.25],'LineWidth',1.2);

    xlabel(ax,'Aligned X position (pixels)');
    ylabel(ax,'Plate ID');
    yticks(ax,uniquePlates);
    title(ax,dataset_label + " — plate-specific X alignment QC", ...
        'Interpreter','none');
    box(ax,'off');
    set(ax,'TickDir','out');

    exportgraphics(fig, ...
        fullfile(output_dir,dataset_label + "_alignment_QC.pdf"), ...
        'ContentType','vector','BackgroundColor','white');
end

%% ========================================================================
% SUMMARY
% ========================================================================

fprintf('\n============================================================\n');
fprintf('XY PREPROCESSING COMPLETE\n');
fprintf('Dataset: %s\n',dataset_label);
fprintf('Independent plates processed: %d\n',height(plate_summary));
fprintf('Tracks retained: %d\n',size(X_data,2));
fprintf('Output directory: %s\n',output_dir);
fprintf('============================================================\n');
