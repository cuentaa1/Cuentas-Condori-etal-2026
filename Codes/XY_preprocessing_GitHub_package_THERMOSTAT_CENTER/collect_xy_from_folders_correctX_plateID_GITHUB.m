function [X_all,Y_all,plate_id,track_id,plate_name, ...
    track_metadata,plate_summary] = ...
    collect_xy_from_folders_correctX_plateID_GITHUB( ...
    folder_list,x_shift_px,min_frames,total_minutes,fps,mat_file_pattern)
% COLLECT_XY_FROM_FOLDERS_CORRECTX_PLATEID_GITHUB
% =========================================================================
% Extract, align, and pool XY trajectories from multiple independent assay
% plates while preserving plate and original track identity.
%
% INPUTS
% ------
% folder_list
%   One folder per independent assay plate.
%
% x_shift_px
%   One translational X offset per plate, in pixels. In the repository-facing
%   preprocessing script this is computed from the measured midpoint between
%   the two fixed thermostats:
%
%       x_shift_px = 1350 - measured_thermostat_center_px
%       X_aligned  = X_raw + x_shift_px
%
%   Positive values therefore move a trajectory toward larger X coordinates.
%   Y coordinates are unchanged.
%
% min_frames
%   Minimum original track length (experiment_1.track(k).npts) required for
%   retention.
%
% total_minutes
%   Length of the common output timeline.
%
% fps
%   Acquisition frame rate.
%
% mat_file_pattern
%   Filename pattern used to identify the cleaned tracking MAT file inside
%   each plate folder. Default: "*postclean_experiment_1*.mat".
%
% OUTPUTS
% -------
% X_all, Y_all
%   nFrames x nTracks matrices. Each column is one retained trajectory.
%   Frames outside the observed lifetime of a track are NaN.
%
% plate_id
%   Numeric plate identifier for every pooled trajectory column.
%   Each input folder is one independent biological plate.
%
% track_id
%   Original track index within experiment_1.track for every output column.
%
% plate_name
%   Original assay-folder name for every output column.
%
% track_metadata
%   One row per retained trajectory, including plate identity, original
%   track index, X shift, original frame range, and original npts.
%
% plate_summary
%   One row per successfully processed assay plate.
%
% RAW TRACKING FORMAT
% -------------------
% Each assay folder must contain a MAT file matching mat_file_pattern and
% containing variable "experiment_1". Coordinates are read from:
%
%       experiment_1.track(k).dq.sloc
%
% where row 1 is X and row 2 is Y. startFrame/endFrame are assumed to use
% the zero-based frame indexing of the tracking output and are converted to
% MATLAB's one-based indexing.
% =========================================================================

    if nargin < 3 || isempty(min_frames)
        min_frames = 300;
    end

    if nargin < 4 || isempty(total_minutes)
        total_minutes = 15;
    end

    if nargin < 5 || isempty(fps)
        fps = 2;
    end

    if nargin < 6 || strlength(string(mat_file_pattern)) == 0
        mat_file_pattern = "*postclean_experiment_1*.mat";
    end

    folder_list = string(folder_list(:));
    x_shift_px = double(x_shift_px(:));

    if numel(folder_list) ~= numel(x_shift_px)
        error('Number of X shifts must match number of plate folders.');
    end

    if any(~isfinite(x_shift_px))
        error('All X shifts must be finite numeric values.');
    end

    if min_frames < 1 || total_minutes <= 0 || fps <= 0
        error('min_frames, total_minutes, and fps must be positive.');
    end

    nFrames = round(total_minutes*60*fps);

    % Trajectory storage.
    Xcols = cell(0,1);
    Ycols = cell(0,1);

    PlateID = zeros(0,1);
    TrackID = zeros(0,1);
    PlateName = strings(0,1);
    XShiftPx = zeros(0,1);
    RawStartFrame = zeros(0,1);
    RawEndFrame = zeros(0,1);
    OriginalNpts = zeros(0,1);
    NStoredFrames = zeros(0,1);

    % Plate-level audit information.
    SummaryPlateID = zeros(0,1);
    SummaryPlateName = strings(0,1);
    SummaryFolder = strings(0,1);
    SummaryMatFile = strings(0,1);
    SummaryXShiftPx = zeros(0,1);
    SummaryNTracksOriginal = zeros(0,1);
    SummaryNTracksRetained = zeros(0,1);

    for fidx = 1:numel(folder_list)

        fldr = folder_list(fidx);
        dx = x_shift_px(fidx);

        if ~isfolder(fldr)
            warning('Plate folder not found; skipping:\n%s',fldr);
            continue
        end

        matInfo = dir(fullfile(fldr,mat_file_pattern));

        if isempty(matInfo)
            warning('No file matching "%s" found in:\n%s', ...
                mat_file_pattern,fldr);
            continue
        end

        % Ambiguous folders should be fixed rather than silently selecting
        % an arbitrary tracking file.
        if numel(matInfo) > 1
            names = string({matInfo.name});
            error(['More than one file matching "%s" was found in:\n%s\n' ...
                'Matches:\n%s\nEach plate folder must resolve to one file.'], ...
                mat_file_pattern,fldr,strjoin(names,newline));
        end

        matPath = fullfile(string(matInfo.folder),string(matInfo.name));

        S = load(matPath,'experiment_1');

        if ~isfield(S,'experiment_1')
            warning('Variable experiment_1 not found; skipping:\n%s',matPath);
            continue
        end

        E = S.experiment_1;

        if ~local_has_member(E,'track')
            warning('experiment_1.track is missing; skipping:\n%s',matPath);
            continue
        end

        T = local_get_member(E,'track');

        if isempty(T)
            warning('experiment_1.track is empty; skipping:\n%s',matPath);
            continue
        end

        [~,folderName] = fileparts(fldr);
        nRetainedThisPlate = 0;

        for k = 1:numel(T)

            Tk = T(k);

            requiredTrackMembers = ["npts","startFrame","endFrame","dq"];
            if ~all(arrayfun(@(q) local_has_member(Tk,q),requiredTrackMembers))
                continue
            end

            npts = double(local_get_member(Tk,'npts'));

            if isempty(npts) || ~isfinite(npts) || npts < min_frames
                continue
            end

            dq = local_get_member(Tk,'dq');

            if ~local_has_member(dq,'sloc')
                continue
            end

            sloc = double(local_get_member(dq,'sloc'));

            if isempty(sloc) || size(sloc,1) < 2
                continue
            end

            startF = double(local_get_member(Tk,'startFrame'));
            endF = double(local_get_member(Tk,'endFrame'));

            if isempty(startF) || isempty(endF) || ...
                    ~isfinite(startF) || ~isfinite(endF) || endF < startF
                continue
            end

            % Raw tracking frames are zero-based. Convert to MATLAB indexing.
            frameVec = (startF:endF) + 1;

            % Coordinate arrays and frame vector are occasionally not
            % identical in length. Retain only the shared available range.
            L = min(size(sloc,2),numel(frameVec));

            if L < 1
                continue
            end

            frameVec = frameVec(1:L);
            xRaw = sloc(1,1:L);
            yRaw = sloc(2,1:L);

            % -------------------------------------------------------------
            % PLATE ALIGNMENT
            % -------------------------------------------------------------
            xAligned = xRaw + dx;

            % Keep only samples that fall inside the requested common assay
            % timeline. Missing times remain NaN.
            withinTimeline = frameVec >= 1 & frameVec <= nFrames;

            if ~any(withinTimeline)
                continue
            end

            frameKeep = frameVec(withinTimeline);

            colX = NaN(nFrames,1);
            colY = NaN(nFrames,1);

            colX(frameKeep) = xAligned(withinTimeline);
            colY(frameKeep) = yRaw(withinTimeline);

            Xcols{end+1,1} = colX; %#ok<AGROW>
            Ycols{end+1,1} = colY; %#ok<AGROW>

            PlateID(end+1,1) = fidx; %#ok<AGROW>
            TrackID(end+1,1) = k; %#ok<AGROW>
            PlateName(end+1,1) = string(folderName); %#ok<AGROW>
            XShiftPx(end+1,1) = dx; %#ok<AGROW>
            RawStartFrame(end+1,1) = startF; %#ok<AGROW>
            RawEndFrame(end+1,1) = endF; %#ok<AGROW>
            OriginalNpts(end+1,1) = npts; %#ok<AGROW>
            NStoredFrames(end+1,1) = numel(frameKeep); %#ok<AGROW>

            nRetainedThisPlate = nRetainedThisPlate + 1;
        end

        SummaryPlateID(end+1,1) = fidx; %#ok<AGROW>
        SummaryPlateName(end+1,1) = string(folderName); %#ok<AGROW>
        SummaryFolder(end+1,1) = fldr; %#ok<AGROW>
        SummaryMatFile(end+1,1) = matPath; %#ok<AGROW>
        SummaryXShiftPx(end+1,1) = dx; %#ok<AGROW>
        SummaryNTracksOriginal(end+1,1) = numel(T); %#ok<AGROW>
        SummaryNTracksRetained(end+1,1) = nRetainedThisPlate; %#ok<AGROW>

        fprintf('Plate %d/%d: %s | shift %+g px | retained %d/%d tracks\n', ...
            fidx,numel(folder_list),folderName,dx, ...
            nRetainedThisPlate,numel(T));
    end

    if isempty(Xcols)
        X_all = NaN(nFrames,0);
        Y_all = NaN(nFrames,0);
    else
        X_all = cat(2,Xcols{:});
        Y_all = cat(2,Ycols{:});
    end

    plate_id = PlateID;
    track_id = TrackID;
    plate_name = PlateName;

    track_metadata = table( ...
        PlateID,PlateName,TrackID,XShiftPx, ...
        RawStartFrame,RawEndFrame,OriginalNpts,NStoredFrames, ...
        'VariableNames',{ ...
        'PlateID','PlateName','OriginalTrackID','XShift_px', ...
        'RawStartFrame_zeroBased','RawEndFrame_zeroBased', ...
        'OriginalNpts','NFramesStored'});

    plate_summary = table( ...
        SummaryPlateID,SummaryPlateName,SummaryFolder,SummaryMatFile, ...
        SummaryXShiftPx,SummaryNTracksOriginal,SummaryNTracksRetained, ...
        'VariableNames',{ ...
        'PlateID','PlateName','Folder','TrackingMATFile','XShift_px', ...
        'NTracksOriginal','NTracksRetained'});

    fprintf('\nPooling complete: %d retained tracks from %d processed plates.\n', ...
        size(X_all,2),height(plate_summary));
end


function tf = local_has_member(value,name)
% True for a struct field or object property.

    name = char(name);

    if isstruct(value)
        tf = isfield(value,name);
    elseif isobject(value)
        tf = isprop(value,name);
    else
        tf = false;
    end
end


function out = local_get_member(value,name)
% Retrieve a struct field or object property.

    name = char(name);

    if isstruct(value)
        out = value.(name);
    elseif isobject(value)
        out = value.(name);
    else
        error('Requested member "%s" from unsupported value type.',name);
    end
end
