function model = process_model_linear(varargin)
% --- process_model() -----------------------------------------------------
% Returns a process model defined by F and Qs. Default is a CV model for a 
% 4D state vector with Ts = 1.
%
% 2023-10-30 Robin Forsling

model = read_inputs(nargin,varargin);

% --- CHOOSE PROCESS MODEL ------------------------------------------------
T = model.Ts; q = model.q; d = model.ncoord;

% CONSTANT POSITION MODEL
if strcmpi(model.model_type,'cp')
    F = eye(d);
    Q = q^2*T*eye(d);

% CONSTANT VELOCITY MODEL
elseif strcmpi(model.model_type,'cv')
    F = [eye(d) T*eye(d) ; zeros(d) eye(d)];
    Q = q^2*[T^3*eye(d)/3 T^2*eye(d)/2 ; T^2*eye(d)/2 T*eye(d)];

% CONSTANT ACCELERATION MODEL
elseif strcmpi(model.model_type,'ca')
    F = [eye(d) T*eye(d) T^2*eye(d)/2 ; zeros(d) eye(d) T*eye(d) ; zeros(d) zeros(d) eye(d)];
    Q = q^2*[T^5*eye(d)/20 T^4*eye(d)/8 T^3*eye(d)/6 ; T^4*eye(d)/8 T^3*eye(d)/3 T^2*eye(d)/2 ; T^3*eye(d)/6 T^2*eye(d)/2 T*eye(d)];
end
% -------------------------------------------------------------------------

model.F = F;
model.Q = Q;

end
% -------------------------------------------------------------------------


% --- READ INPUTS ---------------------------------------------------------
function model = read_inputs(NARG,VARGIN)

    model = process_model_parameters;
  
    model.model_type = 'cv';                                                % Type of model (CV or CA)
    model.ncoord = 2;                                                       % Number of spatial coordinates
    model.nx = 4;                                                           % State dimensionality
    model.Ts = 1;                                                           % Sampling time (optional)
    model.q = 1;
    model.v_max = 100;
    model.a_max = 10;

    if NARG > 0
        if isstruct(VARGIN{1})
            inparams = VARGIN{1};
            if isfield(inparams,'ncoord'); model.ncoord = inparams.ncoord; end
            if isfield(inparams,'Ts'); model.Ts = inparams.Ts; end
            if isfield(inparams,'q'); model.q = inparams.q; end
            if isfield(inparams,'model_type'); model.model_type = inparams.model_type; elseif isfield(inparams,'model'); model.model_type = inparams.model; end
            if isfield(inparams,'v_max'); model.v_max = inparams.v_max; end
            if isfield(inparams,'a_max'); model.a_max = inparams.a_max; end
        else
            if NARG > 0 && ~isempty(VARGIN{1}); model.ncoord = VARGIN{1}; end
            if NARG > 1 && ~isempty(VARGIN{2}); model.Ts = VARGIN{2}; end
            if NARG > 2 && ~isempty(VARGIN{3}); model.model_type = VARGIN{3}; end
        end
    end

    if strcmpi(model.model_type,'cp') || strcmpi(model.model_type,'cpm')
        model.model_type = 'cp';
        if ~isempty(model.ncoord); model.nx = model.ncoord; 
        else; error('invalid dimensionality')
        end
    elseif strcmpi(model.model_type,'cv') || strcmpi(model.model_type,'cvm')
        model.model_type = 'cv';
        if ~isempty(model.ncoord); model.nx = 2*model.ncoord; 
        else; error('invalid dimensionality')
        end
    elseif strcmpi(model.model_type,'ca') || strcmpi(model.model_type,'cam')
        model.model_type = 'ca';
        if ~isempty(model.ncoord); model.nx = 3*model.ncoord; 
        else; error('invalid unknown dimensionality')
        end
    else
        error('unknown model type')
    end

end
% -------------------------------------------------------------------------

