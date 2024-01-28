function par = merge_motion_primitives(varargin)
% --- merge_motion_primitives() -------------------------------------------
% Merges motion primitives defined according to motion_primitive_parameters
%
% Function call:
% merge_motion_primitives(prim1,prim2,prim3,...)
% where primi can be either a motion primitive or a cell array of motion
% primitives.
%
% 2023-10-30 Robin Forsling

par = [];

[motion_prim,nprim] = extract_motion_primitives(nargin,varargin);

if nprim >= 1
    par = motion_prim{1};
    ncoord = size(par.X,1);
    for k = 2:nprim
        temp = motion_prim{k}; N = temp.N;
        if size(temp.X,1) ~= ncoord; error('merge_motion_primitives: incompatible motion primitives...'); end
        dL = norm(par.X(1:2,end)-temp.X(1:2,1));                            % Position difference between last point in current primitive and first point in the next primitive
        if dL <= 1e-3
            K = par.N;
            par.X = [par.X temp.X(:,2:N)]; par.head = [par.head temp.head(:,2:N)]; 
            par.X(3:4,K) = (par.X(3:4,K-1)+par.X(3:4,K+1))/2;
            par.X(5:6,K) = (par.X(5:6,K-1)+par.X(5:6,K+1))/2;
            par.L = par.L + temp.L; par.Ld = par.Ld + temp.Ld;
        else
            par.X = [par.X temp.X]; par.head = [par.head temp.head]; 
            par.L = par.L + temp.L + dL; par.Ld = par.Ld + temp.Ld + dL;
        end
        par.N = size(par.X,2);
    end
end

end


function [motion_prim,nprim] = extract_motion_primitives(NARG,VARGIN)
    
    % Compute number of motion primitives and validate
    nprim = 0;
    for k = 1:NARG
        if iscell(VARGIN{k})
            mpa = VARGIN{k}; nmpa = length(mpa); 
            for l = 1:nmpa
                if validate_motion_primitive_parameters(mpa{l}); nprim = nprim + 1;
                else; error('merge_motion_primitives: invalid motion primitive...')
                end
            end
        else
            if validate_motion_primitive_parameters(VARGIN{k}); nprim = nprim + 1;
            else; error('merge_motion_primitives: invalid motion primitive...')
            end
        end
    end
    if nprim == 0; motion_prim = []; return; end

    % Extract motion primitives
    motion_prim = cell(nprim,1); iprim = 1;
    for k = 1:NARG
        if iscell(VARGIN{k})
            mpa = VARGIN{k}; nmpa = length(mpa); 
            for l = 1:nmpa
                motion_prim{iprim} = mpa{l}; iprim = iprim + 1;
            end
        else
            motion_prim{iprim} = VARGIN{k}; iprim = iprim + 1;
        end
    end
end