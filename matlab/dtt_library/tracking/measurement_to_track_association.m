function A = measurement_to_track_association(association_model,sensor_model,tracks,meas,varargin)
% --- measurement_to_track_association() ----------------------------------
% Measurement to track association.
%
% 2023-10-30 Robin Forsling


% HANDLE INPUTS
if isempty(meas); A = []; return; end
    
% RUN ASSOCIATION
switch association_model.method
    case 0; A = bypass_association(tracks,meas);
%     case 1; A = nn_association(tracks,meas,sensor_model,params);
%     case 2; A = gnn_association(tracks,meas,sensor_model,params);
    otherwise; error('unknown association logic')
end

end
% -------------------------------------------------------------------------


% --- BYPASS ASSSOCIATION -------------------------------------------------
function A = bypass_association(tracks,meas)
    ntracks = length(tracks); nmeas = length(meas);
    A = NaN(ntracks,2); A(:,1) = (1:ntracks)';
    for i = 1:ntracks
        for j = 1:nmeas
            if meas{j}.detection && (meas{j}.ID == tracks{i}.ID)
                A(i,2) = j;
                break;
            end
        end
    end
end
% -------------------------------------------------------------------------


% --- MAHALANOBIS DISTANCE ------------------------------------------------
function d2 = mahalanobis_distance()
    d2 = Inf;
end
% -------------------------------------------------------------------------


% --- LOG-LIKELIHOOD RATIO ------------------------------------------------
function llr = loglikelihood_ratio()
    llr = Inf;
end
% -------------------------------------------------------------------------

