function tracker = create_single_target_tracker(varargin)
% --- create_single_target_tracker() --------------------------------------
% Creates a single sensor tracker able to track a single target.
%
% VARARGIN
%   process_model = varargin{1}
%   sensor_model = varargin{2}
%   fusion_model = varargin{3}
%   association_model = varargin{4}
%   datalink_model = vargargin{5}
%
% 2023-10-30 Robin Forsling

models = read_inputs(nargin,varargin);

tracker = tracker_model_parameters;
tracker.params = tracker_parameters;
tracker.type = 'STT';
tracker.nsensor = 1;
tracker.ntracks = 1;
tracker.process_model = models.process;
tracker.sensor_model = models.sensor;
tracker.fusion_model = models.fusion;
tracker.association_model = models.asso;
tracker.datalink_model = models.datalink;
tracker.tracks = {track_parameters};

tracker.tracks{1}.ID = 1;


end


% --- READ INPUTS ---------------------------------------------------------
function models = read_inputs(NARG,VARGIN)

    % VARARGINS
    if NARG > 0 && ~isempty(VARGIN{1}); process_model = VARGIN{1}; else; process_model = process_model_parameters; end
    if NARG > 1 && ~isempty(VARGIN{2}); sensor_model = VARGIN{2}; else; sensor_model = sensor_model_parameters; end
    if NARG > 2 && ~isempty(VARGIN{3}); fusion_model = VARGIN{3}; else; fusion_model = fusion_model_parameters; end
    if NARG > 3 && ~isempty(VARGIN{4}); association_model = VARGIN{4}; else; association_model = association_model_parameters; end
    if NARG > 4 && ~isempty(VARGIN{5}); datalink_model = VARGIN{5}; else; datalink_model = datalink_model_parameters; end

    % VALIDATE
    if size(sensor_model.pos,1) ~= process_model.ncoord; error('invalid sensor position'); end

    % SET MODELS
    models.process = process_model;
    models.sensor = sensor_model;
    models.fusion = fusion_model;
    models.asso = association_model;
    models.datalink = datalink_model;
end

