function A = track_to_track_association(association_model,fusion_model,tracks1,tracks2,varargin)
% --- track_to_track_association() ----------------------------------------
% Track to track association.
%
% VARARGIN
%   time/params = varargin{1}: current simulated time/tracker params
%
% 2023-10-30 Robin Forsling


% HANDLE INPUTS
if isempty(tracks2); A = []; return; end
if ~iscell(tracks1); tracks1 = {tracks1}; end
if ~iscell(tracks2); tracks2 = {tracks2}; end

% RUN ASSOCIATION
switch association_model.method
    case 0; A = bypass_association(tracks1,tracks2);
%     case 1; A = nn_association(tracks,meas,sensor_model,params);
%     case 2; A = gnn_association(tracks,meas,sensor_model,params);
    otherwise; error('unknown association logic')
end

end


% --- BYPASS ASSSOCIATION -------------------------------------------------
function A = bypass_association(tracks1,tracks2)
    ntracks1 = length(tracks1); ntracks2 = length(tracks2);
    A = NaN(ntracks1,2); A(:,1) = (1:ntracks1)';
    for i = 1:ntracks1
        for j = 1:ntracks2
            if tracks1{i}.init && tracks2{j}.init && (tracks2{j}.ID == tracks1{i}.ID)
                A(i,2) = j;
                break;
            end
        end
    end
end

