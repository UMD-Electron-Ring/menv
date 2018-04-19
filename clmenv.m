% clmenv = command line menv. This is an edit to Hui Li's GUI, which is a
% copy of SPOT by Chris Allen. 
%
% Notes: 
% The GUI is updated to include reference trajectories, this change is not
% reflected in clmenv yet.
% Otherwise, clmenv is much more advanced than Gui menv, changes need to be
% propagated through
%
% Kiersten R. 
% 3/7/17
%

% -- main, initializing function
function clm = clmenv()

clm.open = @openspt;
clm.saveas = @savefile;

clm.defparam = @defparam; % define beam, sim. params
clm.defmatcher = @defmatcher; % load optimization settings
clm.deflattice = @deflattice; % load lattice definition
clm.draw = @drawlattice; % plot lattice structure on figure
clm.build = @buildFODO; % create a basic FODO lattice.

% -- integrate env. for defined problem
clm.solve = @solvelattice;

% -- optimization routines
clm.periodicmatcher = @periodicmatcher; % target is matched periodic solution
clm.targetmatcher = @targetmatcher; % target is x,y,x',y' [optional D,D',nux,nuy,ref-trajectory]
clm.GSmatcher = @GStargetmatcher; % GlobalSearch method, requires globalsearch optimization toolbox
clm.MOmatcher = @MOtargetmatcher; % multi-objective search method, still not functional 


% -- set-up for optimizations
clm.maketarget = @maketarget;
clm.makeoptiset = @makeoptiset;
clm.maketraj = @maketraj; % load ref-trajectory data, plot on figure
clm.makeparams = @makeparams;

clm.thisFig = figure(); axes();
end

% -- functions to be used from command line

function openspt(varargin)
% Open .spt file describing lattice
% --  optional argument is string with filename
global clm

% -- Import file
if nargin == 1
    filename = varargin{1};
    p = genpath(pwd); addpath(p);
    fullpathname = which(filename);
else
    [filename, pathname] = uigetfile('*.spt', 'Open File');
    if( filename==0 ) return; end
    len = length( filename );
    if( len<5 || ~strcmpi(filename(len-3:len),'.spt') )
        errordlg( 'The filename must have an extension .spt!', 'Error', 'modal');
        return;
    end
    fullpathname = [pathname filename];
end


fcopyfile( fullpathname, 'savetmp.mat' );
load 'savetmp';
clm.usrdata = usrdata;
clm.usrdata.filename = fullpathname;  % important to rename its filename


% -- the rest is ensuring backwards compatibility with old .spt files from
% Hui Li / MENV examples folder.
% -- Dipole index:
if( isfield(usrdata,'did')==0 )
    clm.usrdata.did = zeros( size(usrdata.opt) );
end

% -- need to move ic, targets, weights to sub-structures (new format)
if not(exist('clm.usrdata.ic'))
clm.usrdata.ic = struct();
clm.usrdata.ic.x0 = clm.usrdata.x0;
clm.usrdata.ic.xp0 = clm.usrdata.xp0;
clm.usrdata.ic.y0 = clm.usrdata.y0;
clm.usrdata.ic.yp0 = clm.usrdata.yp0;
clm.usrdata=rmfield(clm.usrdata,{'x0','xp0','y0','yp0'});
end

if not(exist('clm.usrdata.target'))
clm.usrdata.target = struct();
clm.usrdata.target.x1 = clm.usrdata.x1;
clm.usrdata.target.xp1 = clm.usrdata.xp1;
clm.usrdata.target.y1 = clm.usrdata.y1;
clm.usrdata.target.yp1 = clm.usrdata.yp1;
clm.usrdata=rmfield(clm.usrdata,{'x1','xp1','y1','yp1'});
end

if not(exist('clm.usrdata.weights'))
clm.usrdata.weights = struct();
clm.usrdata.weights.xw = clm.usrdata.xw;
clm.usrdata.weights.xpw = clm.usrdata.xpw;
clm.usrdata.weights.yw = clm.usrdata.yw;
clm.usrdata.weights.ypw = clm.usrdata.ypw;
clm.usrdata=rmfield(clm.usrdata,{'xw','xpw','yw','ypw'});
end

% -- Transfer from SI
clm.usrdata = TransferFromSI( clm.usrdata );

end

function savefile( filename )
% Save clmenv memory to named file
clobal clm
usrdata = clm.usrdata;
save 'savetmp' usrdata
fcopyfile( 'savetmp.mat', filename );
delete( 'savetmp.mat' );
end

function defparam()
% Load initial beam parameters into clmenv memory
global clm
load 'params'

clm.usrdata.emitance = params.emitance;
clm.usrdata.perveance = params.perveance;
if isfield(params,'ic')
    clm.usrdata.ic = params.ic;
    if not(isfield(clm.usrdata.ic,'D0'))
        clm.usrdata.ic.D0 = 0;
        clm.usrdata.ic.Dp0 = 0;
        warning('Dispersion not specified, setting D=D''=0')
    end
else
    clm.usrdata.ic = struct();
    clm.usrdata.ic.x0 = params.x0;
    clm.usrdata.ic.y0 = params.y0;
    clm.usrdata.ic.xp0 = params.xp0;
    clm.usrdata.ic.yp0 = params.yp0;
    try
        clm.usrdata.ic.D0 = params.D0;
        clm.usrdata.ic.Dp0 = params.Dp0;
    catch
        clm.usrdata.ic.D0 = 0;
        clm.usrdata.ic.Dp0 = 0;
    end
end
clm.usrdata.stepsize = params.stepsize;
clm.usrdata.distance = params.distance;
try clm.usrdata.d0 = params.d0;
catch clm.usrdata.d0 = 0;
end
end

function defmatcher()
global clm

% need to make target file earlier. Target file has structure target.<name>
% which specifies targets and their weights.
% optiset.tmp has settings for optimizer
load 'target'
load 'optiset'

clm.usrdata.maxIter = optiset.maxIter;
clm.usrdata.tolFun = optiset.tolFun;

% -- load targets
clm.usrdata.target = struct();
clm.usrdata.target.x1 = target.x1;
clm.usrdata.target.y1 = target.y1;
clm.usrdata.target.xp1 = target.xp1;
clm.usrdata.target.yp1 = target.yp1;
% -- try adding dispersion targets, set to 0 if undefined
try
    clm.usrdata.target.D1 = target.D1;
    clm.usrdata.target.Dp1 = target.Dp1;
catch
    clm.usrdata.target.D1 = 0;
    clm.usrdata.target.Dp1 = 0;
end

% -- load weights
clm.usrdata.weights = struct();
clm.usrdata.weights.xw = target.xw;
clm.usrdata.weights.yw = target.yw;
clm.usrdata.weights.xpw = target.xpw;
clm.usrdata.weights.ypw = target.ypw;
% -- try adding dispersion weights, set to 0 if undefined
try
    clm.usrdata.weights.Dw = target.Dw;
    clm.usrdata.weights.Dpw = target.Dpw;
catch
    clm.usrdata.weights.Dw = 0;
    clm.usrdata.weights.Dpw = 0;
end

% -- try adding tune targets + weights, set to 0 if undefined
try clm.usrdata.target.nuxt = target.nuxt;
catch clm.usrdata.target.nuxt = 0; end
try clm.usrdata.target.nuyt = target.nuyt;
catch  clm.usrdata.target.nuyt = 0; end
try clm.usrdata.weights.nuxw = target.nuxw;
catch clm.usrdata.weights.nuxw = 0; end
try clm.usrdata.weights.nuyw = target.nuyw;
catch clm.usrdata.weights.nuyw = 0; end

% -- try adding envelope weights, set to 0 if undefined
try clm.usrdata.weights.betaw = target.betaw;
catch clm.usrdata.weights.betaw =0 ; end
try clm.usrdata.weights.refw = target.refw;
catch clm.usrdata.weights.refw =0 ; end


end

function maketarget(varargin)
% Load beam targets and weights into clmenv memory
% varargin:
% 1 -- targetlist, matrix w/ 6 entries, x,y,x',y',D,D' (D,D' optional)
% 2 -- weightlist, matrix w/ 6 entries, x,y,x',y',D,D' (D,D' optional)
% 3 -- optional targetlist, matrix w/ 2 entries nux,nuy
% 4 -- optional weightlist nux,nuy,beta (for symmetry condition betax=betay)

% -- read in target list (4 or 6 entries required)
targetlist = varargin{1};
target.x1 = targetlist(1);
target.y1 = targetlist(2);
target.xp1 = targetlist(3);
target.yp1 = targetlist(4);
try
    target.D1 = targetlist(5);
    target.Dp1 = targetlist(6);
catch
    target.D1 = 0;
    target.Dp1 = 0;
    warning('Dispersion target not specified, set D=D''=0')
end

% -- read in weight list (4 or 6 entries required)
weightlist = varargin{2};
target.xw = weightlist(1);
target.yw = weightlist(2);
target.xpw = weightlist(3);
target.ypw = weightlist(4);
try
    target.Dw = weightlist(5);
    target.Dpw = weightlist(6);
catch
    target.Dw = 0; target.Dpw = 0;
    warning('Dispersion matching will not be included')
end

% -- parse in optional target/weight lists if given. 
if nargin == 2 % no optional targets/weights set
    target.nuxt = 0;
    target.nuyt = 0;
    target.nuxw = 0;
    target.nuyw = 0;
    target.betaw = 0;
    target.refw = 0;
elseif nargin == 3 % envelope weights specified (no tune target/weight)
    optweightlist = varargin{3};
    target.nuxt = 0;
    target.nuyt = 0;
    target.nuxw = 0;
    target.nuyw = 0;
    target.betaw = optweightlist(1);
    target.refw = optweightlist(2);
elseif nargin ==4 % some combination of env + tune target/ weights
    opttargetlist = varargin{3};
    optweightlist = varargin{4};
    if isempty(opttargetlist) && length(optweightlist)==2
        target.nuxt = 0;
        target.nuyt = 0;
        target.nuxw = 0;
        target.nuyw = 0;
        target.betaw = optweightlist(1);
        target.refw = optweightlist(2);
    elseif length(opttargetlist) ==2
        target.nuxt = opttargetlist(1);
        target.nuyt = opttargetlist(2);
        if length(optweightlist)==2
            target.nuxw = optweightlist(1);
            target.nuyw = optweightlist(2);
            target.betaw = 0;
            target.refw = 0;
        elseif length(optweightlist)==3
            target.nuxw = optweightlist(1);
            target.nuyw = optweightlist(2);
            target.betaw = optweightlist(3);
        elseif length(optweightlist)==4
            target.nuxw = optweightlist(1);
            target.nuyw = optweightlist(2);
            target.betaw = optweightlist(3);
            target.refw = optweightlist(4);
        else error('Wrong number of entries in optional list of weights sent to maketarget')
        end
    else error('Wrong number of entries in optional list of targets sent to maketarget')
    end  
end

save 'target' target 
end

function makeoptiset(iterations,tolerance)
% Load iteration and tolerance setting into clmenv memory
optiset.maxIter = iterations;
optiset.tolFun = tolerance;
save 'optiset' optiset 
end

function makeparams(emit,perv,ic,stepsize,startend)
% Make a params file from params arguments

params.emitance = emit; % I know this is misspelled, but if I change it here I have to change it everywhere
params.perveance = perv;
params.stepsize = stepsize;

% -- initial conditions
params.ic = ic;


if length(startend)==1
   params.d0 = 0;
   params.distance = startend;
elseif length(startend)==2
   params.d0 = startend(1);
   params.distance = startend(2);
end

save 'params' params
end

function maketraj(xref,yref)
% Load trajectory data into clmenv memory
global clm

clm.usrdata.target.xref = xref;
clm.usrdata.target.yref = yref;

% -- plot trajectory
if isfield(clm.usrdata,'handle')
if not(isempty(clm.usrdata.handle(1)))
    if ishandle(clm.usrdata.handle(1)) delete(clm.usrdata.handle(1)); end
end
end
hold on
href = plot(xref,yref,':k');
hold off
clm.usrdata.handle(1) = href;

end

function cleartraj()
% Delete stored trajectory data
global clm

% -- delete data
clm.usrdata = rmfield(clm.usrdata,{'xref','yref'});

% -- deleted plotted trajectory
if isfield(clm.usrdata,'handle')
if not(isempty(clm.usrdata.handle(1)))
    if ishandle(clm.usrdata.handle(1)) delete(clm.usrdata.handle(1)); end
end
end

end

function deflattice(varargin)
global clm

if length(varargin)==1
    lat = varargin{1};
    clm.usrdata.ele = lat.ele;
    clm.usrdata.loc = lat.loc;
    clm.usrdata.len = lat.len;
    clm.usrdata.str = lat.str;
    clm.usrdata.opt = lat.opt;
    clm.usrdata.did = lat.did;
    try
        clm.usrdata.irho = lat.irho;
    catch
        clm.usrdata.irho = 0*lat.str;
        warning('Bending radius not defined for elements')
    end
elseif length(varargin)==6
    clm.usrdata.ele = varargin{1};
    clm.usrdata.loc = varargin{2};
    clm.usrdata.len = varargin{3};
    clm.usrdata.str = varargin{4};
    clm.usrdata.opt = varargin{5};
    clm.usrdata.did = varargin{6};
    try
        clm.usrdata.irho = varargin{7};
    catch
        clm.usrdata.irho = 0*varargin{4};
        warning('Bending radius not defined for elements')
    end
end

end

function drawlattice(varargin)
global clm

if isempty(varargin)
    ele = clm.usrdata.ele;
    loc = clm.usrdata.loc;
    len = clm.usrdata.len;
    str = clm.usrdata.str;
    opt = clm.usrdata.opt;
elseif length(varargin)==1
    lat = varargin{1};
    ele = lat.ele;
    loc = lat.loc;
    len = lat.len;
    str = lat.str;
    opt = lat.opt;
elseif length(varargin)==5
    ele = varargin{1};
    loc = varargin{2};
    len = varargin{3};
    str = varargin{4};
    opt = varargin{5};
end

% -- clear figure handle
if isfield(clm.usrdata,'handle')
    if not(isempty(clm.usrdata.handle(2)))
        for i=2:length(clm.usrdata.handle)
            if ishandle(clm.usrdata.handle(i)) delete(clm.usrdata.handle(i)); end
        end
    end
end

% -- plot lattice elements
for i=1:length(ele)
    hele = .25 * str(i)/200;
    %if str(i)==0; continue % -- if str==0, don't draw a patch
    %end
    sele = loc(i)-len(i)*0.5;
    eele = loc(i)+len(i)*0.5;
    % -- choose color based on element, optimization
    if strcmp(ele(i),'D') col = [0,1,0];
    elseif strcmp(ele(i),'Q') && opt(i)==0 col = [.9,.9,.9]; 
    elseif strcmp(ele(i),'Q') && opt(i)==1 col = [.3,.3,.3]; 
    elseif strcmp(ele(i),'S') col = [255,140,0]/255;
    end
    % -- draw a patch for each element
    hpatch = patch([sele,eele,eele,sele],[0,0,hele,hele],col);
    clm.usrdata.handle(i+1) = hpatch;
end
end

function [lat] = buildFODO(ncells)

% -- to-do: make input parser for optional FODO parameters

Iquad = 1.826;
nD = ncells;
nQ = 2*ncells;
nele = nD+nQ;
lcell = 32;
L = lcell*ncells;

% -- dipo params
dlen = 3.819;
dang = 7.8*pi/180;
rho = dlen/dang;
dint = 19.917;
ridg = 338.859;
Idipo = dang*ridg /dint;
dstr = Current2Kappa(Idipo,'D');

% -- quad params
qlen = 5.164; % based on Santiago 2006 Tech. Note
Iquad = 2.194; % Amps, nom. operating pt.
qstr = Current2Kappa(Iquad,'Q');

ele ='';
loc = zeros(1,nele);
len = zeros(1,nele);
str = zeros(1,nele);
did = zeros(1,nele);
opt = zeros(1,nele);
irho = zeros(1,nele);

% -- make list of elements
for i=1:ncells
    ele = [ele,'QDQ'];
end

% -- define quads
for i=1:nQ
    loc(i) = 8 + 16*(i-1);
    len(i) = qlen;
    str(i) = qstr * (-2*mod(i,2)+1);
    did(i) = 0;
    opt(i) = 0;
end

% -- define dipoles
for i=1:nD
    loc(i+nQ) = 16 + 32*(i-1);
    len(i+nQ) = dlen;
    str(i+nQ) = dstr;
    did(i+nQ) = 0.0; % deprecated
    opt(i+nQ) = 0;  
    irho(i+nQ) = 1/rho;
end

% -- sort list
[loc,isort] = sort(loc);

len = len(isort);
str = str(isort);
did = did(isort);
opt = opt(isort);
irho = irho(isort);

lat.ele = ele;
lat.len = len;
lat.loc = loc;
lat.str = str;
lat.did = did;
lat.opt = opt;
lat.irho = irho;

end


function solvelattice()

global clm
runtmp = clm.usrdata;

% check that initial conditions (params) are defined
if( isempty(runtmp.ic.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end

% -- Transfer to SI
runtmp = Transfer2SI( runtmp );

% -- Save to tempary file
save 'runtmp' runtmp;

% -- Run ......
[x,y,xp,yp,D,Dp,d,nux,nuy,Cx,Cy] = runmenv();
% -- transfer results from SI units to cm;
x=x*1.e2; y=y*1.e2; d=d*1.e2;  % m-cm
D=D*1.e2;

% -- Plotting
% -- clear plots if they exist
if isfield(clm,'soldata')
    if isfield(clm.soldata,'handle')
        if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end
        if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end
        if ishandle(clm.soldata.handle(3)) delete(clm.soldata.handle(3)); end
    end
end

hold on; h1 = plot(d,x,'b'); h2 = plot(d,y,'r'); h3 = plot(d,D,'k'); hold off;
xlabel('z (cm)'); ylabel('X:blue, Y:red, D:black (cm)');
axis([ min(d) max(d) 0.0 max([x,y])*1.2 ]);


% Save data, first down-convert resolution if needed.
if( runtmp.stepsize<0.005 ) interval = round(0.005/runtmp.stepsize);
else interval = 1; end
n = length(d); ind = 1:interval:n;
if( ind(length(ind))~=n ) ind(length(ind)+1) = n; end

% -- save plot handles
clm.soldata.handle(1) = h1; 
clm.soldata.handle(2) = h2;
clm.soldata.handle(3) = h3;
% -- save history data 
clm.soldata.d = d(ind); 
clm.soldata.x = x(ind);
clm.soldata.y = y(ind);
clm.soldata.xp = xp(ind);
clm.soldata.yp = yp(ind);
clm.soldata.D = D(ind);
clm.soldata.Dp = Dp(ind);
% -- save tunes + chromaticity
clm.soldata.nux = nux; 
clm.soldata.nuy = nuy;
clm.soldata.Cx = Cx;
clm.soldata.Cy = Cy;
% -- save final conditions
clm.soldata.xf = x(end);
clm.soldata.yf = y(end);
clm.soldata.xpf = xp(end);
clm.soldata.ypf = yp(end);
clm.soldata.Df = D(end);
clm.soldata.Dpf = Dp(end);

end

function periodicmatcher()
% Run matching algorithm for periodic solution for given lattice function

global clm
runtmp = clm.usrdata;
% only check one parameter is enough
if( isempty(runtmp.ic.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end

% Transfer to SI
runtmp = Transfer2SI( runtmp );

% Save to tempary file
save 'runtmp' runtmp;

% Run ......
newX0 = match2period();

% Save the new result
runtmp.ic = newX0;

runtmp = TransferFromSI( runtmp );

clm.usrdata = runtmp;
% Update figure
solvelattice()
end

function targetmatcher()
% Run matching algorithm, find lattice function for desired target (final)
% condition given initial condition
global clm
runtmp = clm.usrdata;
 
% only check one parameter is enough
if( isempty(runtmp.ic.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( isempty(runtmp.target.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( sum(runtmp.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end
   if ishandle(clm.soldata.handle(3)) delete(clm.soldata.handle(3)); end
end

% Transfer to SI
runtmp = Transfer2SI( runtmp );

% Save to temporary file
save 'runtmp' runtmp;

% Run ......
newKappa = match2target;

load runtmp % this might be problematic, moves unwanted data to usrdata?

% Save the new result
[~,n] = size( runtmp.loc ); k = 1;
for i=1:n
    if( runtmp.opt(i) )
        runtmp.str(i) = newKappa(k); k = k+1;
    end
end


runtmp = TransferFromSI( runtmp );
clm.usrdata = runtmp;

% Update figure
solvelattice()
end


function GStargetmatcher()
% Run matching algorithm, vary lattice function to match target, including
% reference trajectory
global clm
runtmp = clm.usrdata;
 
% only check one parameter is enough
if( isempty(runtmp.ic.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( isempty(runtmp.target.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( sum(runtmp.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end
end

% Transfer to SI
runtmp = Transfer2SI( runtmp );

% Save to tempary file
save 'runtmp' runtmp;

% Run ......
[newKappa] = GSmatch2target( 'runtmp' );

% Save the new result
[~,n] = size( runtmp.loc ); k = 1;
for i=1:n
    if( runtmp.opt(i) )
        runtmp.str(i) = newKappa(k); k = k+1;
    end
end

runtmp = TransferFromSI( runtmp );
clm.usrdata = runtmp;

% Update figure
solvelattice()
end

function MOtargetmatcher()
% Run matching algorithm, vary lattice function to match target, including
% reference trajectory
global clm
runtmp = clm.usrdata;
 
% only check one parameter is enough
if( isempty(runtmp.ic.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( isempty(runtmp.target.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end
if( sum(runtmp.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end
end

% Transfer to SI
runtmp = Transfer2SI( runtmp );

% Save to tempary file
save 'runtmp' runtmp;

% Run ......
[newKappa] = MOmatch2target( 'runtmp' );

% Save the new result
[~,n] = size( runtmp.loc ); k = 1;
for i=1:n
    if( runtmp.opt(i) )
        runtmp.str(i) = newKappa(k); k = k+1;
    end
end

runtmp = TransferFromSI( runtmp );
clm.usrdata = runtmp;

% Update figure
solvelattice()
end


% -- sub-functions only ever called from within clmenv()

function usrdata = zeroAxesUserData
usrdata = struct( ...
   'handle',[0,0,0,0,0,0], ...
   'd',[], ...
   'x',[], ...
   'y',[] );
end

function usrdata = zeroMainUserData
usrdata = struct( ...
   'flag',[], ...
   'filename',[], ...
   'emitance',[], ...
   'perveance',[], ...
   'x0',[], ...
   'y0',[], ...
   'xp0',[], ...
   'yp0',[], ...
   'stepsize',0.02, ...
   'distance',[], ...
   'ele',[], ...
   'loc',[], ...
   'len',[], ...
   'str',[], ...
   'did',[], ...
   'opt',[], ...
   'x1',[], ...
   'y1',[], ...
   'xp1',[], ...
   'yp1',[], ...   
   'xw',1.0, ...
   'yw',1.0, ...
   'xpw',1.0, ...
   'ypw',1.0, ...
   'xref',[], ...
   'yref',[], ...
   'maxIter',20, ...
   'tolFun',1.e-8,...
   'handle',[]);
end 

function axesHandle = findAxes( fig )
axesHandle = findobj( fig, 'Type', 'axes' );
for i=1:length(axesHandle)
   pos = get( axesHandle(i), 'Position' );
   x(i) = pos(1);
end;
[x,ind] = sort(x);
% Sort axes handle according to their x position
axesHandle = axesHandle( ind );
end

function usrdata = TransferFromSI( data )
usrdata = data;
usrdata.ic.x0 = data.ic.x0*1.e2;               % m->cm
usrdata.ic.y0 = data.ic.y0*1.e2;               % m->cm
try
    usrdata.ic.D0 = data.ic.D0*1.e2;
catch
    warning('Dispersion D0 not defined, setting D0 = Dp0 = 0')
    usrdata.ic.D0 = 0;
    usrdata.ic.Dp0 = 0;
end
%usrdata.xp0,usrdata.yp0
usrdata.emitance = data.emitance*1.e6;   % mrad->mm.mrad
%usrdata.perveance
usrdata.stepsize = data.stepsize*1.e2;   % m->cm
usrdata.distance = data.distance*1.e2;   % m->cm
try
usrdata.d0 = data.d0*1.e2;   % m->cm
catch
    warning('Start-point s0 not defined, setting s0 = 0')
    usrdata.d0 = 0;
end

%usrdata.ele
usrdata.loc = data.loc*1.e2;             % m->cm
usrdata.len = data.len*1.e2;             % m->cm
%usrdata.str
if isfield(data,'target')
    usrdata.target.x1 = data.target.x1*1.e2;               % m->cm
    usrdata.target.y1 = data.target.y1*1.e2;               % m->cm
end
if isfield(data,'xref')
    usrdata.xref = data.xref*1.e2;               % m->cm
    usrdata.yref = data.yref*1.e2;               % m->cm
end
%usrdata.xp1,usrdata.yp1,usrdata.xw,usrdata.yw,usrdata.xpw,usrdata.ypw,usrdata.maxIter,usrdata.tolFun
end

function usrdata = Transfer2SI( data )
usrdata = data;
usrdata.ic.x0 = data.ic.x0*1.e-2;               % cm->m
usrdata.ic.y0 = data.ic.y0*1.e-2;               % cm->m
try
usrdata.ic.D0 = data.ic.D0*1e-2;
catch end
%usrdata.xp0,usrdata.yp0,usrdata.D0 
usrdata.emitance = data.emitance*1.e-6;   % mm.mrad->mrad
%usrdata.perveance
usrdata.stepsize = data.stepsize*1.e-2;   % cm->m
usrdata.distance = data.distance*1.e-2;   % cm->m
usrdata.d0 = data.d0*1.e-2;   % cm->m
%usrdata.ele
usrdata.loc = data.loc*1.e-2;             % cm->m
usrdata.len = data.len*1.e-2;             % cm->m   
%usrdata.str
if isfield(data,'target')
usrdata.target.x1 = data.target.x1*1.e-2;               % cm->m
usrdata.target.y1 = data.target.y1*1.e-2;               % cm->m
end
if isfield(data,'xref')
usrdata.xref = data.xref*1.e-2;               % m->cm
usrdata.yref = data.yref*1.e-2;               % m->cm
end
%usrdata.xp1,usrdata.yp1,usrdata.xw,usrdata.yw,usrdata.xpw,usrdata.ypw,usrdata.maxIter,usrdata.tolFun
end

function fcopyfile( src, dst )  %force copying src to dst regardless if it exist
% Test if dst exist. If it does, delete it
fid = fopen( dst );
if( fid~=-1 )
   fclose( fid );
   delete( dst );
end
% After delete dst, src can safely overwrite dst.
% Otherwise in UNIX, copyfile(src,dst) will stop until Enter key is pressed in the Command Window.
copyfile( src, dst );
end
