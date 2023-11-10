function tracker = tracker_model_parameters()
% --- tracker_model_parameters() ------------------------------------------
% Creates parameter struct for a generic tracker.
%
% 2023-10-30 Robin Forsling

% DESCRIPTIONS
tracker.type = 'type of tracker (optional), e.g.: STT, MTT, MSMTT';
tracker.nsensor = 'number of sensors, e.g.: 1 for a STT and MTT';
tracker.ntracks = 'number of tracks';
tracker.process_model = 'process model, same is used for all tracks';
tracker.sensor_model = 'sensor model, can differ';
tracker.fusion_model = 'parameters for the fusion functionality';
tracker.association_model = 'parameters for the association functionality';
tracker.datalink_model = 'parameters for the datalink functionality';
tracker.tracks = 'single track or array of tracks';

% SET ALL FIELDS EMPTY
tracker = set_struct_fields_empty(tracker);