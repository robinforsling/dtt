function varargout = sensor_model_spherical(varargin)
% --- sensor_model_spherical() --------------------------------------------
% Returns a sensor model with observations given in polar coordinates 
% (default) or spherical coordinates (m = 3).
%
% 2023-10-30 Robin Forsling

[model,struct_req] = read_inputs(nargin,varargin);
if struct_req; varargout{1} = model; return; end

if nargout == 1
    varargout{1} = model;
else
    varargout{1} = model.h;
    varargout{2} = model.J;
    varargout{3} = model.R;
    varargout{4} = model.PD;
    varargout{5} = model.beta_FA;
end

end
% -------------------------------------------------------------------------


% --- READ INPUTS ---------------------------------------------------------
function [model,struct_req] = read_inputs(NARG,VARGIN)

    [h,J] = transform_cart_to_pol;
    sr = 50; saz = d2r(1);    

    model = sensor_model_parameters;

    model.sensor_type = 'polar';
    model.ori = 0;
    model.fov = 2*pi;
    model.r_lim = [1 1e9];
    model.m = 2;
    model.h = h;
    model.J = J;
    model.R = diag([sr saz].^2);
    model.PD = 1.0;
    model.beta_FA = 0;

    struct_req = 0;
    if NARG == 0
        return;
    elseif NARG == 1
        if func_params_struct_req(VARGIN{1})
            model = get_sensor_model_params('descr');
            struct_req = 1;
            return;
        end
    end

    model.R = [];

    if isstruct(VARGIN)
        inparams = VARGIN{1};
        if isfield(inparams,'m'); model.m = inparams.m; end
        %if isfield(inparams,'R'); model.R = inparams.R; end
        if isfield(inparams,'sigma_r'); sr = inparams.sigma_r; elseif isfield(inparams,'sr'); sr = inparams.sr; end
        if isfield(inparams,'sigma_az'); saz = inparams.sigma_az; elseif isfield(inparams,'saz'); saz = inparams.saz; 
        elseif isfield(inparams,'sigma_ang'); saz = inparams.sigma_ang; elseif isfield(inparams,'sang'); saz = inparams.sang; 
        elseif isfield(inparams,'sigma_phi'); saz = inparams.sigma_phi; elseif isfield(inparams,'sphi'); saz = inparams.sphi; 
        end
        if isfield(inparams,'PD'); model.PD = inparams.PD; end
        if isfield(inparams,'beta_FA'); model.beta_FA = inparams.beta_FA; end
    else
        if NARG > 0 && ~isempty(VARGIN{1}); model.m = VARGIN{1}; end
        if NARG > 1 && ~isempty(VARGIN{2}); sr = VARGIN{2}; end
        if NARG > 2 && ~isempty(VARGIN{3}); saz = VARGIN{3}; end
        if NARG > 3 && ~isempty(VARGIN{4}); model.PD = VARGIN{4}; end
        if NARG > 4 && ~isempty(VARGIN{5}); model.beta_FA = VARGIN{5}; end
    end

    if model.m == 3
        [h,J] = transform_cart_to_sph;
        model.sensor_type = 'spherical';
        model.h = h;
        model.J = J;
    end

    if model.m ~= 2 && model.m ~= 3; error('m must be either 2 or 3...'); end

    if isempty(model.R)
        model.R = diag([sr saz*ones(1,model.m-1)].^2);
    end
end
% -------------------------------------------------------------------------

