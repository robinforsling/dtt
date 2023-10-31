function tracks = fuse_tracks(fusion_model,tracks_own,tracks_rx,params,varargin)
% --- fuse_tracks() -------------------------------------------------------
% Fuse tracks. Track list tracks_dl is allowed to be empty by which it
% follows that the fusing is trivial. Track list tracks_own empty yields a
% warning. For fusion to work, tracks_own and tracks_rx must be of equal 
% length. For multitrack fusion an assignment matrix A must be provided.
%
% VARARGIN
%   A = varargin{1}:
%   
% 2023-10-30 Robin Forsling


% HANDLE INPUTS
if nargin > 5; A = varargin{1}; end
if iscell(tracks_own); is_cell = 1; else; is_cell = 0; tracks_own = {tracks_own}; end
if ~iscell(tracks_rx); tracks_rx = {tracks_rx}; end


% CHECK TRACK LISTS
if isempty(tracks_rx) || isempty(tracks_rx{1}); tracks = tracks_own; return;
elseif isempty(tracks_own); tracks = tracks_own; warning('tracks_own is empty'); return;
elseif length(tracks_own) ~= length(tracks_rx); error('input track lists must be of equal length')
end
if length(tracks_own) > 1 && nargin < 6; error('must provide an assignment matrix A'); end
if length(tracks_own) == 1; A = [1 1]; end 


% FUSE TRACKS
if fusion_model.method == 0                                                 % No fusion
    tracks = tracks_own;                                                    
else                                                                        % Fusion
    ntracks = length(tracks_own); tracks = tracks_own;
    for i = 1:ntracks
        if ~isnan(A(i,2))
            track = tracks_own{i}; 
            track_rx = tracks_rx{A(i,2)}; 
            track_com = track.common_est;
            if track.init && track_rx.dl_est.valid               
                [xho,Po,Ho] = extract_estimate(track);                      % Own estimate
                [xhdl,Pdl,Hdl] = extract_estimate(track_rx.dl_est);         % DL estimates
                [xhc,Pc,Hc] = extract_estimate(track_com);                  % Common estimate
                switch fusion_model.method
                    case 1 % KF
                        [xhat,P] = kf_fusion(xho,Po,Ho,xhdl,Pdl,Hdl);
                        [ghat,G] = kf_fusion(xhc,Pc,Hc,xhdl,Pdl,Hdl); 
                    case 2 % GIMF, uses naive/kf
                        [xhat,P] = kf_fusion(xho,Po,Ho,xhdl,Pdl,Hdl);
                        [ghat,G] = kf_fusion(xhc,Pc,Hc,xhdl,Pdl,Hdl); 
                    case 3 % BSC
                        % Not implemented yet
                    case 4 % CI
                        [xhat,P,w] = ci_fusion(xho,Po,Ho,xhdl,Pdl,Hdl);
                        [ghat,G] = ci_fusion(xhc,Pc,Hc,xhdl,Pdl,Hdl,w);
                    case 5 % ICI
                        % Not implemented yet
                    case 6 % LE
                        [xhat,P,Ko,Kdl,P12] = le_fusion(xho,Po,xhdl,Pdl,Hdl);
                        [ghat,G] = le_fusion(xhc,Pc,xhdl,Pdl,Hdl,Ko,Kdl,P12);
                    case 7 % CI DCA-HYP
                        [xhat,P,w_vec] = ci_dcahyp_fusion(xho,Po,xhdl,Pdl,Hdl);
                        [ghat,G] = ci_dcahyp_fusion(xhc,Pc,xhdl,Pdl,Hdl,w_vec);
                    otherwise; error('unknown fusion method')
                end
                tracks{i} = update_track_and_common_est(track,xhat,P,ghat,G,params);
            elseif ~track.init && track_rx.init                              % Initialize track from tracks_rx
                tracks{i} = initialize_track_using_dl_track(track,track_rx);
            end
        end
    end
end

if ~is_cell; tracks = tracks{1}; end

end


% --- UPDATE TRACK --------------------------------------------------------
function track = update_track_and_common_est(track,xhat,P,ghat,G,params)
    track.xhat = xhat;
    track.P = P;
    track.time.last_fuse = params.time;
    track.common_est.xhat = ghat;
    track.common_est.P = G;
    track.common_est.time.last_fuse = params.time;
end


% --- INITIALIZE TRACK FROM DL TRACK --------------------------------------
function track = initialize_track_using_dl_track(track,dl_track,params)
    track.init = 1;
    track.xhat = dl_track.xhat;
    track.P = dl_track.P;
    track.time.init = params.time;
    track.common_est.xhat = ghat;
    track.common_est.P = G;
end
