function tracks = generate_datalink_tracks(cntrl,dl_model,fus_mod,tracks,params)
% --- generate_datalink_tracks() ------------------------------------------
% Generate datalink estimates for communication.
%
% 2023-10-30 Robin Forsling

if isempty(cntrl) || isempty(tracks); return; end
if ~cntrl.dl.req_estimates; return; end
if iscell(tracks); is_cell = 1; else; is_cell = 0; tracks = {tracks}; end

m = dl_model.comm_mgmt.ndirections;
ntracks = length(tracks);

for i = 1:ntracks
    Psi = []; yPsi = []; RPsi = [];
    if tracks{i}.init && dl_model.comm_mgmt.method ~= 0
        [y2,R2] = extract_estimate(tracks{i});
        [yg,Rg] = extract_estimate(tracks{i}.common_est);
        if cntrl.cheat.globally_known_est 
            [~,R1] = extract_estimate(cntrl.cheat.tracks{i});             
        else 
            R1 = Rg; %R1s = Rg;
        end
        if dl_model.comm_mgmt.loss_function == 2; [y2s,R2s] = subtract_estimate(y2,R2,yg,Rg); end    

        switch dl_model.comm_mgmt.method
            case 0 % No communication (LKF)
                
            case 1 % Full estimate
                Psi = eye(length(y2));
                if fus_mod.method == 2; y2 = y2s; R2 = R2s; % GIMF
                end

            case 2 % PCO
                Psi = dimred_pco(R2,m);

            case 3 % PARO
                if dl_model.comm_mgmt.loss_function == 1; Psi = dimred_paro(R1,R2,m);
                elseif dl_model.comm_mgmt.loss_function == 2; Psi = dimred_paro(R1,R2s,m);
                elseif dl_model.comm_mgmt.loss_function == 4; Psi = dimred_paro(R1,R2,m,'ci');
                else; warning('unknown dimred loss function')
                end

            case 4 % GEVO
                if dl_model.comm_mgmt.loss_function == 1; Psi = dimred_gevo(R1,R2,m);
                elseif dl_model.comm_mgmt.loss_function == 2; Psi = dimred_gevo(R1,R2s,m,'kf');
                elseif dl_model.comm_mgmt.loss_function == 4; Psi = dimred_gevo(R1,R2,m,'ci');
                elseif dl_model.comm_mgmt.loss_function == 6; Psi = dimred_gevo(R1,R2,m,'le');
                else; warning('unknown dimred loss function')
                end

            case 5 % DCA-EIG
                D2 = dca_eig(R2);

            case 6 % DCA-OPT
                D2 = dca_opt(R2);

            case 7 % DCA-DOM
                D2 = dca_dom(R2);

            case 8 % DCA-HYP
                D2 = diag(diag(R2));

            otherwise; warning('unknown selection method')
        end
    end

    % --- DIMENSION-REDUCTION ---
    if is_element_in_vector(dl_model.comm_mgmt.method,1:4)
        if dl_model.comm_mgmt.loss_function == 2; y2 = y2s; R2 = R2s; end
        
        % UPDATE DL ESTIMATE
        yPsi = Psi*y2; RPsi = make_symmetric(Psi*R2*Psi');
        tracks{i} = update_dl_est(tracks{i},yPsi,RPsi,Psi,params);

        % UPDATE COMMON ESTIMATE
        [ghat,G,Hg] = extract_estimate(tracks{i}.common_est);
        switch fus_mod.method
            case 1; [ghat,G] = kf_fusion(ghat,G,Hg,yPsi,RPsi,Psi);
            case 2; [ghat,G] = kf_fusion(ghat,G,Hg,yPsi,RPsi,Psi);          % GIMF uses KF since common estimate has been subtracted
            case 3; % BSC, not implemented yet
            case 4; [ghat,G] = ci_fusion(ghat,G,Hg,yPsi,RPsi,Psi);
            case 5; % ICI, not implemented yet
            case 6; [ghat,G] = le_fusion(ghat,G,yPsi,RPsi,Psi);  
        end
        tracks{i}.common_est.xhat = ghat;
        tracks{i}.common_est.P = make_symmetric(G);
    end

    % --- DIAGONAL COVARIANCE APPROXIMATION ---
    if is_element_in_vector(dl_model.comm_mgmt.method,5:8)

        H = eye(length(y2));

        % UPDATE DL ESTIMATE
        tracks{i} = update_dl_est(tracks{i},y2,D2,H,params);

        % UPDATE COMMON ESTIMATE
        [ghat,G,Hg] = extract_estimate(tracks{i}.common_est);
        switch fus_mod.method
            case 1; [ghat,G] = kf_fusion(ghat,G,Hg,y2,D2,H);
            case 2; [ghat,G] = kf_fusion(ghat,G,Hg,y2,D2,H);                % GIMF uses KF since common estimate has been subtracted
            case 3; % BSC, not implemented yet
            case 4; [ghat,G] = ci_fusion(ghat,G,Hg,y2,D2,H);
            case 5; % ICI, not implemented yet
            case 6; [ghat,G] = le_fusion(ghat,G,y2,D2,H);
            %case 7; [ghat,G] = ci_dcahyp_fusion(ghat,G,y2,D2,H);
        end
        tracks{i}.common_est.xhat = ghat;
        tracks{i}.common_est.P = make_symmetric(G);
    end
end

if ~is_cell; tracks = tracks{1}; end

end
% -------------------------------------------------------------------------


% --- UPDATE DL ESTIMATE --------------------------------------------------
function track = update_dl_est(track,yH,RH,H,params)
    track.dl_est.valid = 1;
    track.dl_est.xhat = yH;
    track.dl_est.P = RH;
    track.dl_est.H = H;
    track.dl_est.time = params.time;
end
% -------------------------------------------------------------------------

