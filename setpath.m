% -- add path
menvlocation = '.';
menvpath = genpath(menvlocation); addpath(menvpath);

% -- make sure windows aren't docked (menv complains):
set(0,'DefaultFigureWindowStyle' , 'normal')