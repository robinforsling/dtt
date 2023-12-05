function clr = get_paper_colors
% --- get_paper_colors() -------------------------------------------------
% Color definitions used in the paper.
%
% 2023-12-05 Robin Forsling

f = 255;


% --- BASE COLORS ---
clr.red = hex2rgb('D73027');
clr.darkred = hex2rgb('A50026');
clr.orange = hex2rgb('FC8D59');
clr.lightorange = hex2rgb('FDAE61');
clr.yellow = hex2rgb('FEE08B');
clr.lightyellow = hex2rgb('FFFFBF');
clr.darkyellow = clr_strength(clr.yellow,2);
clr.green = hex2rgb('91CF60');
clr.lightgreen = hex2rgb('D9EF8B');
clr.darkgreen = hex2rgb('1A9850');
clr.cyan = hex2rgb('5AB4AC');
clr.lightcyan = hex2rgb('C7EAE5');
clr.darkcyan = hex2rgb('01665E');
clr.blue = hex2rgb('67A9CF');
clr.lightblue = hex2rgb('D1E5F0');
clr.darkblue = hex2rgb('2166AC');
clr.purple = hex2rgb('762A83');
clr.lightpurple = hex2rgb('AF8DC3');
clr.brown = hex2rgb('8C510A');
clr.lightbrown = hex2rgb('D8B365');
clr.beige = hex2rgb('F6E8C3');


% --- AGENTS ---
clr.agent1 = clr.darkyellow;
clr.agent2 = clr.blue;
clr.agent3 = clr.orange;
clr.agent4 = clr.green;
clr.agent5 = clr.lightbrown;
clr.agent6 = clr.cyan;
clr.agent7 = clr.lightpurple;
clr.agent8 = clr.red;


% --- METHODS --
clr.kf = hex2rgb('F16913');
clr.ci = hex2rgb('78C679');
clr.ici = clr.darkyellow;
clr.le = clr_strength(clr.blue,1.5);
%clr.bluestimator = clr.blue;
clr.clue = clr.darkyellow;
clr.bclue = hex2rgb('238443');
clr.lower = clr.blue;
clr.upper = clr.darkyellow;
clr.lkf = 0.25*[1 1 1];
clr.dkf = clr.darkyellow;
clr.naive = clr.red;
clr.ro = clr.cyan;

% DR AQ:
clr.ls = clr.cyan;
clr.sm = clr.orange;
clr.lg = clr.green;
clr.full = 0.15*[1 1 1];
clr.fus = clr.orange;
clr.asso = clr.cyan;

clr.m = cell(4,1);
clr.m{1} = clr.red;
clr.m{2} = clr.darkyellow;
clr.m{3} = clr.darkgreen;
clr.m{4} = clr.darkblue;


% --- SPECTRAL: 7 COLORS ---
clr.spectral7 = { ...  
    [213,62,79]/f ; ...
    [252,141,89]/f ; ...
    [254,224,139]/f ; ...
    [255,255,191]/f ; ...
    [230,245,152]/f ; ...
    [153,213,148]/f ; ...
    [50,136,189]/f};


% --- GRAY SCALE ---
clr.white = [1 1 1];
clr.llllightgray = 0.95*[1 1 1];
clr.lllightgray = 0.90*[1 1 1];
clr.llightgray = 0.85*[1 1 1];
clr.lightgray = 0.75*[1 1 1];
clr.gray = 0.5*[1 1 1];
clr.darkgray = 0.25*[1 1 1];
clr.ddarkgray = 0.15*[1 1 1];
clr.dddarkgray = 0.1*[1 1 1];
clr.ddddarkgray = 0.05*[1 1 1];
clr.black = [0 0 0];

end


function rgb = hex2rgb(h)
    r = hex2dec(h(1:2));
    g = hex2dec(h(3:4));
    b = hex2dec(h(5:6));
    rgb = [r,g,b]/255;
end

function rgb_new = clr_strength(rgb,p)
    e = [1,1,1];
    rgb_new = (1-p)*e + p*rgb;
end