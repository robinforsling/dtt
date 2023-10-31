function SEED = get_seed_from_time()
% --- get_seed_from_time() ------------------------------------------------
% Computes seed based on current time.
%
% 2023-10-30 Robin Forsling

c1 = 1e7; c2 = 1e4;
SEED = round(c1*mod(datenum(datetime),1));
SEED = round(mod(SEED,c2));