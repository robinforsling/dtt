function S = basic_communication_schedule(ntx,N,varargin)
% --- basic_communication_schedule() --------------------------------------
% Returns the basic communication schedule.
%
% 2023-10-30 Robin Forsling

if nargin > 2; offset = varargin{1}; else; offset = 0; end

B = eye(ntx); S = [];
for i = 1:ceil((N+offset)/ntx); S = [S B]; end
S = S(:,1+offset:N+offset);