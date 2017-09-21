% clmenv = command line menv. This is an edit to Hui Li's GUI, which is a
% copy of SPOT by Chris Allen. 
%
% Notes: 
% The GUI is updated to include reference trajectories, this change is not
% reflected in clmenv yet.
% clmenv can import .spt files but cannot define new lattices -- .spt files
% can be created in GUI interface and then loaded by clmenv. A "define
% element" function should be added in the future.
%
% Kiersten R. 
% 3/7/17
%

% -- main
function clm = clmenv()

clm.open = @openspt;
clm.saveas = @savefile;

clm.defparam = @defparam;
clm.defmatcher = @defmatcher;
clm.deflattice = @deflattice;
clm.draw = @drawlattice;
clm.build = @buildFODO;

clm.solve = @solvelattice;
clm.periodicmatcher = @periodicmatcher;
clm.targetmatcher = @targetmatcher; % target is x,y,xp,yp
clm.targetmatcher2 = @targetmatcher2; % target is x,y,xp,yp,nux,nuy,betax~betay
clm.trajmatcher = @trajmatcher; % target x,y,xp,yp plus tries to match trajectory

clm.maketarget = @maketarget;
clm.makeoptiset = @makeoptiset;
clm.maketraj = @maketraj;
clm.makeparams = @makeparams;

clm.thisFig = figure(); axes();
end


% -- functions that can be used externally

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
% To be compatible with previous version file with no usrdata.did
if( isfield(usrdata,'did')==0 )
    clm.usrdata.did = zeros( size(usrdata.opt) );
end
% Transfer from SI
clm.usrdata = TransferFromSI( usrdata );

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
clm.usrdata.x0 = params.x0;
clm.usrdata.y0 = params.y0;
clm.usrdata.xp0 = params.xp0;
clm.usrdata.yp0 = params.yp0;
clm.usrdata.stepsize = params.stepsize;
clm.usrdata.distance = params.distance;
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

clm.usrdata.x1 = target.x1;
clm.usrdata.y1 = target.y1;
clm.usrdata.xp1 = target.xp1;
clm.usrdata.yp1 = target.yp1;
clm.usrdata.xw = target.xw;
clm.usrdata.yw = target.yw;
clm.usrdata.xpw = target.xpw;
clm.usrdata.ypw = target.ypw;

try clm.usrdata.nuxw = target.nuxw; catch; end
try clm.usrdata.nuyw = target.nuyw; catch; end
try clm.usrdata.betaw = target.betaw; catch; end


end

function maketarget(varargin)
% Load beam targets and weights into clmenv memory
% varargin:
% 1 -- targetlist
% 2 -- weightlist
% 3 -- optional weightlist
%
if nargin > 0
    targetlist = varargin{1};
    target.x1 = targetlist(1);
    target.y1 = targetlist(2);
    target.xp1 = targetlist(3);
    target.yp1 = targetlist(4);

    weightlist = varargin{2};
    target.xw = weightlist(1);
    target.yw = weightlist(2);
    target.xpw = weightlist(3);
    target.ypw = weightlist(4);
end
if nargin > 2
    opttargetlist = varargin{3};
    optweightlist = varargin{4};
    target.nuxt = opttargetlist(1);
    target.nuyt = opttargetlist(2);
    target.nuxw = optweightlist(1);
    target.nuyw = optweightlist(2);
    target.betaw = optweightlist(3);
end

save 'target' target 
end

function maketraj(xref,yref)
% Load trajectory data into clmenv memory
global clm

clm.usrdata.xref = xref;
clm.usrdata.yref = yref;

% -- plot trajectory
if exist('clm.usrdata.handle')
    if ishandle(clm.usrdata.handle) delete(clm.usrdata.handle); end
end
hold on
href = plot(xref,yref,':k');
hold off
clm.usrdata.handle = href;

end

function makeoptiset(iterations,tolerance)
% Load iteration and tolerance setting into clmenv memory
optiset.maxIter = iterations;
optiset.tolFun = tolerance;
save 'optiset' optiset 
end

function makeparams(emit,perv,x0,y0,xp0,yp0,stepsize,distance)
% Make a params file from params arguments

params.emitance = emit;
params.perveance = perv;
params.x0 = x0;
params.y0 = y0;
params.xp0 = xp0;
params.yp0 = yp0;
params.stepsize = stepsize;
params.distance = distance;

save 'params' params
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
elseif length(varargin)==6
    clm.usrdata.ele = varargin(1);
    clm.usrdata.loc = varargin(2);
    clm.usrdata.len = varargin(3);
    clm.usrdata.str = varargin(4);
    clm.usrdata.opt = varargin(5);
    clm.usrdata.did = varargin(6);
end

end

function drawlattice(varargin)

if length(varargin)==1
    lat = varargin{1};
    ele = lat.ele;
    loc = lat.loc;
    len = lat.len;
    str = lat.str;
    opt = lat.opt;
elseif length(varargin)==5
    ele = varargin(1);
    loc = varargin(2);
    len = varargin(3);
    str = varargin(4);
    opt = varargin(5);
end



for i=1:length(ele)
    hele = 1;
    if str(i)==0; continue % -- if str==0, don't draw a patch
    end
    sele = loc(i)-len(i)*0.5;
    eele = loc(i)+len(i)*0.5;
    % -- choose color based on element, optimization
    if strcmp(ele(i),'D') col = [0,1,0];
    elseif strcmp(ele(i),'Q') && opt(i)==0 col = [.9,.9,.9]; 
    elseif strcmp(ele(i),'Q') && opt(i)==1 col = [.3,.3,.3]; 
    end
    % -- draw a patch for each element
    patch([sele,eele,eele,sele],[0,0,hele,hele],col)
end
end

function [lat] = buildFODO(ncells)

nD = ncells;
nQ = 2*ncells;
nele = nD+nQ;
lcell = 32;
L = lcell*ncells;


ele ='';
loc = zeros(1,nele);
len = zeros(1,nele);
str = zeros(1,nele);
did = zeros(1,nele);
opt = zeros(1,nele);

% -- make list of elements
for i=1:ncells
    ele = [ele,'QDQ'];
end

% -- define quads
for i=1:nQ
    loc(i) = 8 + 16*(i-1);
    len(i) = 3.6400;
    str(i) = 140 * (-2*mod(i,2)+1);
    did(i) = 0;
    opt(i) = 0;
end

% -- define dipoles
for i=1:nD
    loc(i+nQ) = 16 + 32*(i-1);
    len(i+nQ) = 3.8500;
    str(i+nQ) = 15.7;
    did(i+nQ) = 0.72;
    opt(i+nQ) = 0;  
end

% -- sort list
[loc,isort] = sort(loc);

len = len(isort);
str = str(isort);
did = did(isort);
opt = opt(isort);

lat.ele = ele;
lat.len = len;
lat.loc = loc;
lat.str = str;
lat.did = did;
lat.opt = opt;

end


function solvelattice()

global clm
usrdata = clm.usrdata;

% check that initial conditions (params) are defined
if( isempty(usrdata.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end

% Transfer to SI
usrdata = Transfer2SI( usrdata );

% Save to tempary file
save 'runtmp' usrdata;

% Run ......
[x,y,xp,yp,d,nux,nuy] = runmenv( 'runtmp' );
x=x*1.e2; y=y*1.e2; d=d*1.e2;  % m-cm

% Plot
% -- clear plots if they exist
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end
end

hold on; h1 = plot(d,x,'b'); h2 = plot(d,y,'r'); hold off;
xlabel('z (cm)'); ylabel('X:blue, Y:red (cm)');
axis([ min(d) max(d) 0.0 max([x,y])*1.2 ]);


% Save data
if( usrdata.stepsize<0.005 ) interval = round(0.005/usrdata.stepsize);
else interval = 1; end
n = length(d); ind = 1:interval:n;
if( ind(length(ind))~=n ) ind(length(ind)+1) = n; end

% -- save plot handles
clm.soldata.handle(1) = h1; 
clm.soldata.handle(2) = h2;
% -- save history data 
clm.soldata.d = d(ind); 
clm.soldata.x = x(ind);
clm.soldata.y = y(ind);
% -- save tunes
clm.soldata.nux = nux; 
clm.soldata.nuy = nuy;
% -- save final conditions
clm.soldata.xf = x(end);
clm.soldata.yf = y(end);
clm.soldata.xpf = xp(end);
clm.soldata.ypf = yp(end);

end

function periodicmatcher()
% Run matching algorithm for periodic solution for given lattice function

global clm
usrdata = clm.usrdata;
% only check one parameter is enough
if( isempty(usrdata.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( isempty(usrdata.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end;

% Transfer to SI
usrdata = Transfer2SI( usrdata );

% Save to tempary file
save 'runtmp' usrdata;

% Run ......
newX0 = match2period( 'runtmp' );

% Save the new result
usrdata.x0 = newX0(1);
usrdata.y0 = newX0(2);
usrdata.xp0= newX0(3);
usrdata.yp0= newX0(4);
usrdata = TransferFromSI( usrdata );

clm.usrdata = usrdata;
% Update figure
solvelattice()
end

function targetmatcher()
% Run matching algorithm, find lattice function for desired target (final)
% condition given initial condition
global clm
usrdata = clm.usrdata;
 
% only check one parameter is enough
if( isempty(usrdata.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( isempty(usrdata.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( sum(usrdata.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end;

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end;
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end;
end

% Transfer to SI
usrdata = Transfer2SI( usrdata );

% Save to tempary file
save 'runtmp' usrdata;

% Run ......
newKappa = match2target( 'runtmp' );

% Save the new result
[~,n] = size( usrdata.loc ); k = 1;
for i=1:n
    if( usrdata.opt(i) )
        usrdata.str(i) = newKappa(k); k = k+1;
    end;
end;

usrdata = TransferFromSI( usrdata );
clm.usrdata = usrdata;

% Update figure
solvelattice()
end

function targetmatcher2()
% Run matching algorithm, find lattice function for desired target (final)
% condition given initial condition
global clm
usrdata = clm.usrdata;
 
% only check one parameter is enough
if( isempty(usrdata.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( isempty(usrdata.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( sum(usrdata.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end;

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end;
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end;
end

% Transfer to SI
usrdata = Transfer2SI( usrdata );

% Save to tempary file
save 'runtmp' usrdata;

% Run ......
newKappa = match2target2( 'runtmp' );

% Save the new result
[~,n] = size( usrdata.loc ); k = 1;
for i=1:n
    if( usrdata.opt(i) )
        usrdata.str(i) = newKappa(k); k = k+1;
    end;
end;

usrdata = TransferFromSI( usrdata );
clm.usrdata = usrdata;

% Update figure
solvelattice()
end

function trajmatcher()
% Run matching algorithm, vary lattice function to match target, including
% reference trajectory
global clm
usrdata = clm.usrdata;
 
% only check one parameter is enough
if( isempty(usrdata.x0) )
    warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( isempty(usrdata.x1) )
    warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
    return;
end;
if( sum(usrdata.opt)==0 )
    warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
    return;
end;

% clear plot if necessary
if exist('clm.soldata')
   if ishandle(clm.soldata.handle(1)) delete(clm.soldata.handle(1)); end;
   if ishandle(clm.soldata.handle(2)) delete(clm.soldata.handle(2)); end;
end

% Transfer to SI
usrdata = Transfer2SI( usrdata );

% Save to tempary file
save 'runtmp' usrdata;

% Run ......
newKappa = match2reftraj( 'runtmp' );

% Save the new result
[~,n] = size( usrdata.loc ); k = 1;
for i=1:n
    if( usrdata.opt(i) )
        usrdata.str(i) = newKappa(k); k = k+1;
    end;
end;

usrdata = TransferFromSI( usrdata );
clm.usrdata = usrdata;

% Update figure
solvelattice()
end


% -- functions only used with clmenv()

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
usrdata.x0 = data.x0*1.e2;               % m->cm
usrdata.y0 = data.y0*1.e2;               % m->cm
%usrdata.xp0,usrdata.yp0
usrdata.emitance = data.emitance*1.e6;   % mrad->mm.mrad
%usrdata.perveance
usrdata.stepsize = data.stepsize*1.e2;   % m->cm
usrdata.distance = data.distance*1.e2;   % m->cm
%usrdata.ele
usrdata.loc = data.loc*1.e2;             % m->cm
usrdata.len = data.len*1.e2;             % m->cm
%usrdata.str
usrdata.x1 = data.x1*1.e2;               % m->cm
usrdata.y1 = data.y1*1.e2;               % m->cm
if exist('data.xref')
    usrdata.xref = data.xref*1.e2;               % m->cm
    usrdata.yref = data.yref*1.e2;               % m->cm
end
%usrdata.xp1,usrdata.yp1,usrdata.xw,usrdata.yw,usrdata.xpw,usrdata.ypw,usrdata.maxIter,usrdata.tolFun
end

function usrdata = Transfer2SI( data )
usrdata = data;
usrdata.x0 = data.x0*1.e-2;               % cm->m
usrdata.y0 = data.y0*1.e-2;               % cm->m
%usrdata.xp0,usrdata.yp0   
usrdata.emitance = data.emitance*1.e-6;   % mm.mrad->mrad
%usrdata.perveance
usrdata.stepsize = data.stepsize*1.e-2;   % cm->m
usrdata.distance = data.distance*1.e-2;   % cm->m
%usrdata.ele
usrdata.loc = data.loc*1.e-2;             % cm->m
usrdata.len = data.len*1.e-2;             % cm->m   
%usrdata.str
usrdata.x1 = data.x1*1.e-2;               % cm->m
usrdata.y1 = data.y1*1.e-2;               % cm->m
if exist('data.xref')
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