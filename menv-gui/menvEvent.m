function menvEvent( event )
switch( get(gcbf,'Tag') )
case 'MainFig'
   MainFigWndProc( gcbf, event );
case 'DefElementFig'
   DefElementWndProc( gcbf, event );
case 'DefParamFig'
   DefParamWndProc( gcbf, event );
case 'DefMatcherFig'
   DefMatcherWndProc( gcbf, event );
end

function MainFigWndProc( thisFig, event )
global guim
switch( event )
case 'Create'
   usrdata = zeroMainUserData;
   guim = clmenv('nofig'); 
%   set( thisFig, 'UserData', usrdata );
   prepareMenu( thisFig );
case 'Axes1Create'
   axdata = zeroAxesUserData;
   set( gcbo, 'UserData', axdata );
case 'Exit'
   close( thisFig );
case 'About'
   msgbox('MENV-(Matlab Envelope Solver), Author: Hui Li + Kiersten Ruisard', 'ABOUT MENV', 'none', 'modal' );
case 'New'
   menvEvent( 'Close' );
case 'Open'
   guim.open();
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, guim.usrdata, 1 );
   [~,n] = size( guim.usrdata.loc );
   changeButtonState( thisFig, n );
case 'Close'
   % Clear current axis
   axesHandle = findAxes( thisFig );
   clearAxes( axesHandle(1) );
   % Clear userData
   guim.usrdata = zeroMainUserData;
   set( thisFig, 'UserData', usrdata );
   % Clear listbox
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   set( listboxHandle , 'String', '' );
   changeButtonState( thisFig, 0 );
case 'Save'
   %usrdata = get( thisFig, 'UserData' );
   if( exist('guim.userdata.filename') )
      if ( isempty('guim.usrdata.filename') )
          menvEvent( 'SaveAs' );
      end
   else
      % Transfer to SI
      usrdata = Transfer2SI( guim.usrdata );
      savefile( usrdata, guim.usrdata.filename );
   end
case 'SaveAs'
   [filename, pathname] = uiputfile('*.spt', 'Save As');
   if( filename==0 ) return; end
   dot = length( find(filename=='.') );
   if( dot==0 )
      filename=[filename '.spt'];
   else
      len = length( filename );
      if( len<5 || ~strcmpi(filename(len-3:len),'.spt') )
         errordlg( 'The filename must have an extension .spt!', 'Error', 'modal');
         return;
      end
   end
   usrdata = get( thisFig, 'UserData' );
   filename = [pathname filename];
   guim.usrdata.filename = filename;
   set( thisFig, 'UserData', usrdata );
   % Transfer to SI
   usrdata = Transfer2SI( guim.usrdata );   
   savefile( guim.usrdata, filename );
case 'Edit'
   usrdata = get( thisFig, 'UserData' );
   % Pass the calling instruction   
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   item = get( listboxHandle, 'Value' );
   guim.usrdata.flag = item;
   set( thisFig, 'UserData', usrdata );
   % Call the input dialogbox
   fig = defElement;
   set( fig, 'UserData', thisFig );
   % Fill out the input dialogbox
   elePopupmenuHandle = findobj( fig, 'Tag', 'ElementPopupMenu' );
   if( guim.usrdata.ele(item)=='Q' )
      set( elePopupmenuHandle, 'Value', 1 );
   elseif( guim.usrdata.ele(item)=='S' )
      set( elePopupmenuHandle, 'Value', 2 );
   else
      set( elePopupmenuHandle, 'Value', 3 );
   end
   locEditHandle = findobj( fig, 'Tag', 'LocationEditText' );
   lenEditHandle = findobj( fig, 'Tag', 'LengthEditText' );
   strEditHandle = findobj( fig, 'Tag', 'StrengthEditText' );
   irhoEditHandle = findobj( fig, 'Tag', 'IrhoEditText' );
   locString = num2str( guim.usrdata.loc(item) );
   lenString = num2str( guim.usrdata.len(item) );
   strString = num2str( guim.usrdata.str(item) );
   try irhoString = num2str( guim.usrdata.irho(item) ); catch irhoString='0'; end
   set( locEditHandle, 'String', locString );
   set( lenEditHandle, 'String', lenString ); 
   set( strEditHandle, 'String', strString );
   set( irhoEditHandle, 'String', irhoString );
   optCheckboxHandle = findobj( fig, 'Tag', 'OptimCheckbox' );
   set( optCheckboxHandle, 'Value', guim.usrdata.opt(item) );
   % important in here to disable/enable some controls
   DefElementWndProc( fig, 'ElementChange' ); 
case 'Insert'
   usrdata = get( thisFig, 'UserData' );
   guim.usrdata.flag = 0;
   set( thisFig, 'UserData', usrdata );
   fig = defElement;
   set( fig, 'UserData', thisFig );
   % important in here to disable/enable some controls
   DefElementWndProc( fig, 'ElementChange' ); 
case 'Delete'
   usrdata = get( thisFig, 'UserData' );
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   item = get( listboxHandle, 'Value' );
   [~,n] = size( usrdata.loc );
   if( n==1 )
      ind = [];
   else
      tmp = 1:n; tmp(item) = n+1; [~,ind] = sort(tmp);
      ind = ind(1:n-1);
   end
   usrdata.ele = usrdata.ele(ind);
   usrdata.loc = usrdata.loc(ind);
   usrdata.len = usrdata.len(ind);
   usrdata.str = usrdata.str(ind);
   usrdata.irho = usrdata.irho(ind);
   usrdata.opt = usrdata.opt(ind);
   set( thisFig, 'UserData', usrdata );
   % Update the listbox
   updateListboxString( listboxHandle, guim.usrdata, 1 );
   changeButtonState( thisFig, n-1 );
case 'Param'
   fig = defParam;
   %set( fig, 'UserData', thisFig );
   %usrdata = get( thisFig, 'UserData' );
   % Initial conditions
   if( isfield(guim.usrdata,'ic') )
      x0EditHandle = findobj( fig, 'Tag', 'X0EditText' );
      y0EditHandle = findobj( fig, 'Tag', 'Y0EditText' );
      xp0EditHandle = findobj( fig, 'Tag', 'XP0EditText' );
      yp0EditHandle = findobj( fig, 'Tag', 'YP0EditText' );
      D0EditHandle = findobj( fig, 'Tag', 'D0EditText' );
      Dp0EditHandle = findobj( fig, 'Tag', 'DP0EditText' );
      set( x0EditHandle, 'String', num2str(guim.usrdata.ic.x0) );
      set( y0EditHandle, 'String', num2str(guim.usrdata.ic.y0) );
      set( xp0EditHandle, 'String', num2str(guim.usrdata.ic.xp0) );
      set( yp0EditHandle, 'String', num2str(guim.usrdata.ic.yp0) );
      set( D0EditHandle, 'String', num2str(guim.usrdata.ic.D0) );
      set( Dp0EditHandle, 'String', num2str(guim.usrdata.ic.Dp0) );
   end
   % Global beam parameters
   if( isfield(guim.usrdata,'emitance') )
      emitEditHandle = findobj( fig, 'Tag', 'EmitEditText' );
      set( emitEditHandle, 'String', num2str(guim.usrdata.emitance) );
   end
   if( isfield(guim.usrdata,'perveance') )
      pvnEditHandle = findobj( fig, 'Tag', 'PvnEditText' );
      set( pvnEditHandle, 'String', num2str(guim.usrdata.perveance) );
   end
   % Numerical parameters
   if( isfield(guim.usrdata,'stepsize') )
      stepEditHandle = findobj( fig, 'Tag', 'StepEditText' );
      set( stepEditHandle, 'String', num2str(guim.usrdata.stepsize) );
   end
   if( isfield(guim.usrdata,'distance') )
      distEditHandle = findobj( fig, 'Tag', 'DistEditText' );
      set( distEditHandle, 'String', num2str(guim.usrdata.distance) );
   end
case 'Solution'
    guim.solve();
%    runtmp = guim.usrdata; %get( thisFig, 'UserData' );
%    % only check one parameter is enough
%    if( isempty(runtmp.ic.x0) )
%       warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    % Transfer to SI
%    runtmp = Transfer2SI( runtmp );
%    % Save to tempary file
%    save 'runtmp' runtmp;
%    % Run ......
%    oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
%    [x,y,xp,yp,D,Dp,d,nux,nuy,Cx,Cy] = runmenv( 'runtmp' );
%    x=x*1.e2; y=y*1.e2; d=d*1.e2;  % m-cm
%    D=D*1.e2;
%    set( thisFig, oldptr{:} );
%    % Axes
%    axesHandle = findAxes( thisFig );
%    axes( axesHandle(1) );
%    axdata = get( axesHandle(1), 'UserData' );
%    % Plot
%    if( axdata.handle(1)~=0 ) delete(axdata.handle(1)); end
%    if( axdata.handle(2)~=0 ) delete(axdata.handle(2)); end
%    hold on; h1 = plot(d,x,'b'); h2 = plot(d,y,'r'); hold off;
%    xlabel('z (cm)'); ylabel('X:blue, Y:red (cm)');
%    axis([ min(d) max(d) 0.0 max([x,y])*1.2 ]);
%    
%     % Save data
%    if( usrdata.stepsize<0.005 ) interval = round(0.005/usrdata.stepsize);
%    else interval = 1; end
%    n = length(d); ind = 1:interval:n;
%    if( ind(length(ind))~=n ) ind(length(ind)+1) = n; end   
%    axdata.handle(1)=h1; axdata.handle(2)=h2; 
%    axdata.d=d(ind); axdata.x=x(ind); axdata.y=y(ind);
%    axdata.xp=xp(ind); axdata.yp=yp(ind);
%    axdata.nux = nux; axdata.nuy = nuy;
%    set( axesHandle, 'UserData', axdata );   
case 'Draw'
    guim.draw();
    yl = ylim();
    ylim([-.3,yl(2)]);
case 'MatcherParam'
    fig = defMatcher;
    %   set( fig, 'UserData', thisFig );
    %   usrdata = get( thisFig, 'UserData' );
    % Target conditions
    if( isfield(guim.usrdata,'target') )
        if( isfield(guim.usrdata.target,'x1') )
            x1EditHandle = findobj( fig, 'Tag', 'X1EditText' );
            y1EditHandle = findobj( fig, 'Tag', 'Y1EditText' );
            xp1EditHandle = findobj( fig, 'Tag', 'XP1EditText' );
            yp1EditHandle = findobj( fig, 'Tag', 'YP1EditText' );
            nux1EditHandle = findobj( fig, 'Tag', 'nux1EditText' );
            nuy1EditHandle = findobj( fig, 'Tag', 'nuy1EditText' );
            set( x1EditHandle, 'String', num2str(guim.usrdata.target.x1) );
            set( y1EditHandle, 'String', num2str(guim.usrdata.target.y1) );
            set( xp1EditHandle, 'String', num2str(guim.usrdata.target.xp1) );
            set( yp1EditHandle, 'String', num2str(guim.usrdata.target.yp1) );
            try set( nux1EditHandle, 'String', num2str(guim.usrdata.target.nux1) );
            catch set( nux1EditHandle, 'String', '0'); end
            try set( nuy1EditHandle, 'String', num2str(guim.usrdata.target.nuy1) );
            catch set( nuy1EditHandle, 'String', '0'); end
        end
    end
    % Weights
    if( isfield(guim.usrdata,'weights') )
        if( isfield(guim.usrdata.weights,'xw') )
            xwEditHandle = findobj( fig, 'Tag', 'XWEditText' );
            ywEditHandle = findobj( fig, 'Tag', 'YWEditText' );
            xpwEditHandle = findobj( fig, 'Tag', 'XPWEditText' );
            ypwEditHandle = findobj( fig, 'Tag', 'YPWEditText' );
            DwEditHandle = findobj( fig, 'Tag', 'DwEditText' );
            DpwEditHandle = findobj( fig, 'Tag', 'DpwEditText' );
            nuxwEditHandle = findobj( fig, 'Tag', 'nuxwEditText' );
            nuywEditHandle = findobj( fig, 'Tag', 'nuywEditText' );
            betawEditHandle = findobj( fig, 'Tag', 'betawEditText' );
            refwEditHandle = findobj( fig, 'Tag', 'refwEditText' );
            set( xwEditHandle, 'String', num2str(guim.usrdata.weights.xw) );
            set( ywEditHandle, 'String', num2str(guim.usrdata.weights.yw) );
            set( xpwEditHandle, 'String', num2str(guim.usrdata.weights.xpw) );
            set( ypwEditHandle, 'String', num2str(guim.usrdata.weights.ypw) );
            % -- optional weights set to 0 if un-used
            try set( DwEditHandle, 'String', num2str(guim.usrdata.weights.Dw) );
            catch set( DwEditHandle, 'String', 0); end
            try set( DpwEditHandle, 'String', num2str(guim.usrdata.weights.Dpw) );
            catch set( DpwEditHandle, 'String', 0); end
            try set( nuxwEditHandle, 'String', num2str(guim.usrdata.weights.nuxw) );
            catch set( nuxwEditHandle, 'String', 0); end
            try set( nuywEditHandle, 'String', num2str(guim.usrdata.weights.nuyw) );
            catch set( nuywEditHandle, 'String', 0); end
            try set( betawEditHandle, 'String', num2str(guim.usrdata.weights.betaw) );
            catch set( betawEditHandle, 'String', 0); end
            try set( refwEditHandle, 'String', num2str(guim.usrdata.weights.refw) );
            catch set( refwEditHandle, 'String', 0); end
        end
    end
   % iterations
   if( isfield(guim.usrdata,'maxIter') )
      maxIterEditHandle = findobj( fig, 'Tag', 'MaxIterEditText' );
      tolFunEditHandle = findobj( fig, 'Tag', 'TolFunEditText' );
      set( maxIterEditHandle, 'String', num2str(guim.usrdata.maxIter) );
      set( tolFunEditHandle, 'String', num2str(guim.usrdata.tolFun) );
   end
case 'Run Periodic Matcher'
    guim.periodicmatcher();
% %   usrdata = get( thisFig, 'UserData' );
%    % only check one parameter is enough
%    if( isempty(guim.usrdata.ic.x0) )
%       warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    if( isempty(guim.usrdata.target.x1) )
%       warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    % Transfer to SI
%    usrdata = Transfer2SI( guim.usrdata );
%    % Save to tempary file
%    save 'runtmp' usrdata;
%    % Run ......
%    %oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
%    oldptr = get(thisFig,'pointer');  set( thisFig, 'pointer', 'watch' );
%    newX0 = match2period();
%    set( thisFig, 'pointer', oldptr );
%    % Save the new result
%    guim.usrdata.ic.x0 = newX0(1);
%    guim.usrdata.ic.y0 = newX0(2);
%    guim.usrdata.ic.xp0= newX0(3);
%    guim.usrdata.ic.yp0= newX0(4);
%    guim.usrdata = TransferFromSI( guim.usrdata );
%    %   set( thisFig, 'UserData', usrdata );
%    % Update figure
%    menvEvent( 'Solution' );
case 'Matcher'
    guim.targetmatcher();
%    usrdata = get( thisFig, 'UserData' );
%    % only check one parameter is enough
%    if( isempty(guim.usrdata.ic.x0) )
%       warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    if( isempty(guim.usrdata.target.x1) )
%       warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    if( sum(guim.usrdata.opt)==0 )
%       warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
%       return;
%    end
%    % Transfer to SI
%    usrdata = Transfer2SI( guim.usrdata );
%    % Save to tempary file
%    save 'runtmp' usrdata;
%    % Run ......
%    oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
%    newKappa = match2target( 'runtmp' );
%    set( thisFig, oldptr{:} );
%    % Save the new result
%    %usrdata = get( thisFig, 'UserData' );
%    [~,n] = size( usrdata.loc ); k = 1;
%    for i=1:n
%       if( guim.usrdata.opt(i) )
%          guim.usrdata.str(i) = newKappa(k); k = k+1;
%       end
%    end
%    %set( thisFig, 'UserData', usrdata );
   % Update the listbox
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, guim.usrdata, 1 );   
   
case 'Global Search'
   guim.GSmatcher();   
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, guim.usrdata, 1 );  
   
case 'Multi-Objective'
   guim.MOmatcher();
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, guim.usrdata, 1 );  
   
    
case 'Coordinate'
   coordTextHandle = findobj( thisFig, 'Tag', 'CoordStaticText' );
   axesHandle = findAxes( thisFig );
   xlim = get( axesHandle, 'XLim' );
   yl = get( axesHandle, 'YLim' );
   set( coordTextHandle, 'Visible', 'on' );
   while( 1 )
      [x,y] = ginput(1);
      str = sprintf('(%f, %f)',x,y);
      set( coordTextHandle, 'String', str );
      if( x<xlim(1) || x>xlim(2) || y<yl(1) || y>yl(2) ) break; end
   end
   set( coordTextHandle, 'String', '' );
   set( coordTextHandle, 'Visible', 'off' );
case 'Import'
   [filename, pathname] = uigetfile('*.txt', 'Open File');
   if( filename==0 ) return; end
   filename = [pathname filename];
   fid = fopen( filename );
   A = [];
   while( ~feof(fid) )
      str = fgets( fid ); data = sscanf( str, '%f' );
      try
         if( ~isempty(data) && (length(data)==2 || length(data)==3) )
            A = [A data];
         end
      catch
         errordlg( 'File format error while reading!', 'ERROR', 'modal' );
         fclose(fid); return;
      end
   end
   fclose( fid );
   if( isempty(A) )
      errordlg( 'No data to be read!', 'ERROR', 'modal' );
   end
   % Preview data and set options
   fig = ImportOpt;
   listboxHandle = findobj( fig, 'Tag', 'PreviewListbox' );
   [~,n] = size(A);
   listboxString = cell(n,1);
   for i=1:n
      itemString = sprintf('%f  ',A(:,i) );
      listboxString(i) = { itemString };
   end
   set( listboxHandle, 'String', listboxString );
   uiwait( fig );
   if( ishandle(fig) )
      radiobuttonHandle = findobj( fig, 'Style', 'radiobutton', 'Value', 1);
      switch( get(radiobuttonHandle, 'Tag') )
      case 'mmRadiobutton'
         A = A*0.1;
      case 'mRadiobutton'
         A = A*1.e2;
      case 'cmRadiobutton'
      end
      delete(fig);
   else
      % == cm   
   end
   % Plot
   axesHandle = findAxes( thisFig );
   axdata = get( axesHandle(1), 'UserData' );
   axes( axesHandle(1) ); hold on;
   %style = { 'k*', 'mo' };
   style = { 'k.', 'm.' };
   [n,~] = size( A );
   x = A(1,:);
   for i=2:n
      y = A(i,:);
      if( axdata.handle(i+1)~=0 ) delete(axdata.handle(i+1)); end
      h = plot( x, y, style{i-1} );
      axdata.handle(i+1) = h;
   end
   hold off;
   % Save
   set( axesHandle, 'UserData', axdata );
case 'ClearImport'
   axesHandle = findAxes( thisFig );
   axdata = get( axesHandle(1), 'UserData' );
   for i=3:4
      if( axdata.handle(i)~=0 )
         delete( axdata.handle(i) );
         axdata.handle(i) = 0;
      end
   end
   set( axesHandle, 'UserData', axdata );   
case 'Export'
   axesHandle = findAxes( thisFig );
   axdata = get( axesHandle(1), 'UserData' );
   if( ~isempty(axdata) && axdata.handle(1)~=0 && axdata.handle(2)~=0 )
      A = [axdata.d; axdata.x; axdata.y];
   else
      warndlg( 'No figure to be exported!', 'WARNING', 'modal' );
      return;
   end
   [filename, pathname] = uiputfile('*.*', 'Save As');
   if( filename==0 ) return; end
   filename = [pathname filename];
   fid = fopen( filename, 'w' );
   fprintf( fid, 'Dist\tX\tY\r\n' );
   fprintf( fid, '%12.8e\t%12.8e\t%12.8e\r\n',A );
   fclose( fid );
case 'Kappa to Current'
%   usrdata = get( thisFig, 'UserData' );
   % only check one parameter is enough
   if( isempty(guim.usrdata.ic.x0) )
      warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   % Transfer to SI
   usrdata = Transfer2SI( guim.usrdata );
   [~,n] = size( usrdata.loc );
   if( n==0 )
      return;
   end
   title = 'Kappa to Current';
   for i=1:n
      if( usrdata.ele(i)=='Q' )
         I = Kappa2Current( usrdata.str(i) );
         msg{i} = ['Quad: ' num2str(I) 'A'];
      elseif( usrdata.ele(i)=='S' )
         I = Kappa2Current( usrdata.str(i), 'S' );
         msg{i} = ['Sol:  ' num2str(I) 'A'];
      else
         I = 0.0;
         msg{i} = ['Dipl:  ' 'Unknown' 'A'];
      end
   end
   listdlg( 'PromptString','Element Current:', ...
            'SelectionMode','single',...
            'ListSize', [160 160], ...
            'ListString',msg);
        
case 'Draw Reference Trajectory'
    %usrdata = get( thisFig, 'UserData' );
    % Draw on figure
    axesHandle = findAxes( thisFig );
    [xref,yref] = DrawReferenceTrajectory(axesHandle);
    % save in UserData
    guim.usrdata.target.xref = xref;
    guim.usrdata.target.yref = yref;
    % set( thisFig, 'UserData', usrdata );
    % Save new reference trajectory to file
    [filename, pathname] = uiputfile('*.txt', 'Save As');
    dlmwrite([pathname,filename],[xref,yref]);
    % Plot on axes
    hold on
    if isfield(guim.usrdata,'href'), delete(guim.usrdata.href); end
    guim.usrdata.href = plot(xref,yref,':k');
    
case 'Load Reference Trajectory'
    %usrdata = get( thisFig, 'UserData' );    
    % load .txt file
    [filename, pathname] = uigetfile('*.txt', 'Open File');
    if( filename==0 ) return; end
    len = length( filename );
    filename = [pathname filename];
    temp =  importdata(filename);
    xref = temp(:,1); yref = temp(:,2);    
    % save in UserData
    guim.usrdata.target.xref = xref;
    guim.usrdata.target.yref = yref;
    %set( thisFig, 'UserData', usrdata );
    % Plot on axes
    hold on
    if isfield(guim.usrdata,'href'), delete(guim.usrdata.href); end
    guim.usrdata.href = plot(xref,yref,':k');

otherwise
end


function DefElementWndProc( thisFig, event )
global guim
switch( event )
case 'OK'
   % Judge if all the inputs valid
   elePopupmenuHandle = findobj( thisFig, 'Tag', 'ElementPopupMenu' );
   locEditHandle = findobj( thisFig, 'Tag', 'LocationEditText' );
   lenEditHandle = findobj( thisFig, 'Tag', 'LengthEditText' );
   strEditHandle = findobj( thisFig, 'Tag', 'StrengthEditText' );
   irhoEditHandle = findobj( thisFig, 'Tag', 'IrhoEditText' );
   element = get( elePopupmenuHandle, 'Value' );
   locString = get( locEditHandle, 'String' );
   lenString = get( lenEditHandle, 'String' );
   strString = get( strEditHandle, 'String' );
   irhoString = get( irhoEditHandle, 'String' );
   location = str2double( locString );
   length = str2double( lenString );
   strength = str2double( strString );
   invrho = str2double( irhoString );
   if( isnan(location) || isnan(length) || isnan(strength) || (element==3 && isnan(invrho)) )
      errordlg( 'Some input parameters are not numbers!', 'ERROR', 'modal' );
      return;
   end
   optCheckboxHandle = findobj( thisFig, 'Tag', 'OptimCheckbox' );
   optim = get( optCheckboxHandle, 'Value' );
   % Find the shared userdata
   mainFigHandle = get( thisFig, 'UserData' );
   listboxHandle = findobj( mainFigHandle, 'Tag', 'ElementListbox' );
%   usrdata = get( mainFigHandle, 'UserData' );
   % 1:Quad; 2:Sol; 3:Dipl
   if( element==1 )
      element = 'Q';
      invrho = 0;
   elseif( element==2 )
      element = 'S';
      invrho = 0;
   else
      element = 'D';
      optim = 0; % disable optim for Dipl
   end
   % 0:Insert; 1:Edit
   if( guim.usrdata.flag==0 )
      guim.usrdata.ele = [guim.usrdata.ele element];
      guim.usrdata.loc = [guim.usrdata.loc location];
      guim.usrdata.len = [guim.usrdata.len length];
      guim.usrdata.str = [guim.usrdata.str strength];
      guim.usrdata.irho = [guim.usrdata.irho invrho];
      guim.usrdata.opt = [guim.usrdata.opt optim];
   else
      item = guim.usrdata.flag;
      guim.usrdata.ele(item) = element;
      guim.usrdata.loc(item) = location;
      guim.usrdata.len(item) = length;
      guim.usrdata.str(item) = strength;
      guim.usrdata.irho(item) = invrho;
      guim.usrdata.opt(item) = optim;
   end
   % Find who is the item that we are editing or inserting
   [~,n] = size( guim.usrdata.loc ); 
   [tmp,ind] = sort( guim.usrdata.loc );
   if( guim.usrdata.flag==0 )
      focus = find( ind==n );
   else
      focus = find( ind==item );
   end
   % Sort
   guim.usrdata.loc = tmp; 
   guim.usrdata.ele = guim.usrdata.ele(ind);
   guim.usrdata.len = guim.usrdata.len(ind);
   guim.usrdata.str = guim.usrdata.str(ind);
   guim.usrdata.irho = guim.usrdata.irho(ind);
   guim.usrdata.opt = guim.usrdata.opt(ind);
   % Save
%   set( mainFigHandle, 'UserData', usrdata );
   % Add formated string to the caller listbox
   updateListboxString( listboxHandle, guim.usrdata, focus );
   changeButtonState( mainFigHandle, n );
   close( thisFig );
case 'Cancel'
   close( thisFig );
case 'ElementChange'
   elePopupmenuHandle = findobj( thisFig, 'Tag', 'ElementPopupMenu' );
   irhoEditHandle = findobj( thisFig, 'Tag', 'IrhoEditText' );
   optCheckboxHandle = findobj( thisFig, 'Tag', 'OptimCheckbox' );
   element = get( elePopupmenuHandle, 'Value' );
   if( element==3 ) % Dipl
      set( optCheckboxHandle, 'Enable', 'off' );
      set( irhoEditHandle, 'Enable', 'on' );
      set( irhoEditHandle, 'BackgroundColor', [1 1 1] ); % white
   else
      set( optCheckboxHandle, 'Enable', 'on' );
      set( irhoEditHandle, 'Enable', 'off' );
      set( irhoEditHandle, 'BackgroundColor', [0.7529 0.7529 0.7529] ); % gray      
   end
end


function DefParamWndProc( thisFig, event )
global guim
switch( event )
case 'OK'
   % Initial conditions
   x0EditHandle = findobj( thisFig, 'Tag', 'X0EditText' );
   y0EditHandle = findobj( thisFig, 'Tag', 'Y0EditText' );
   xp0EditHandle = findobj( thisFig, 'Tag', 'XP0EditText' );
   yp0EditHandle = findobj( thisFig, 'Tag', 'YP0EditText' ); 
   D0EditHandle = findobj( thisFig, 'Tag', 'D0EditText' );
   Dp0EditHandle = findobj( thisFig, 'Tag', 'DP0EditText' ); 
   x0 = str2double( get( x0EditHandle, 'String' ) );
   y0 = str2double( get( y0EditHandle, 'String' ) );
   xp0 = str2double( get( xp0EditHandle, 'String' ) );
   yp0 = str2double( get( yp0EditHandle, 'String' ) );   
   D0 = str2double( get( D0EditHandle, 'String' ) );
   Dp0 = str2double( get( Dp0EditHandle, 'String' ) );   
   if( isnan(x0) || isnan(y0) || isnan(xp0) || isnan(yp0) || isnan(D0) || isnan(Dp0) ||x0<=0 || y0<=0 )
      errordlg( 'Initial condition input error!', 'ERROR', 'modal' );
      return;   
   end
   
   % Global beam parameters
   emitEditHandle = findobj( thisFig, 'Tag', 'EmitEditText' );
   pvnEditHandle = findobj( thisFig, 'Tag', 'PvnEditText' );
   emitance = str2double( get( emitEditHandle, 'String' ) );
   perveance = str2double( get( pvnEditHandle, 'String' ) );
   if( isnan(emitance) || isnan(perveance ) || emitance<=0 )
      errordlg( 'Global beam parameters input error!', 'ERROR', 'modal' );
      return;
   end
   % Numerical parameters
   stepEditHandle = findobj( thisFig, 'Tag', 'StepEditText' );
   distEditHandle = findobj( thisFig, 'Tag', 'DistEditText' );
   stepsize = str2double( get( stepEditHandle, 'String' ) );
   distance = str2double( get( distEditHandle, 'String' ) );   
   if( isnan(stepsize) || isnan(distance) || stepsize<=0 || distance<=0 )
      errordlg( 'Numerical parameters input error!', 'ERROR', 'modal' );
      return;
   end
   % Find the shared userdata
%   mainFigHandle = get( thisFig, 'UserData' );
%   usrdata = get( mainFigHandle, 'UserData' );
   guim.usrdata.emitance = emitance;
   guim.usrdata.perveance = perveance;
   guim.usrdata.ic.x0 = x0;
   guim.usrdata.ic.y0 = y0;
   guim.usrdata.ic.xp0 = xp0;
   guim.usrdata.ic.yp0 = yp0;
   guim.usrdata.ic.D0 = D0;
   guim.usrdata.ic.Dp0 = Dp0;
   guim.usrdata.stepsize = stepsize;
   guim.usrdata.distance = distance;
%   set( mainFigHandle, 'UserData', usrdata );
case 'Cancel'
end
close( thisFig );

function DefMatcherWndProc( thisFig, event )
global guim
switch( event )
case 'OK'
   % Target conditions
   x1EditHandle = findobj( thisFig, 'Tag', 'X1EditText' );
   y1EditHandle = findobj( thisFig, 'Tag', 'Y1EditText' );
   xp1EditHandle = findobj( thisFig, 'Tag', 'XP1EditText' );
   yp1EditHandle = findobj( thisFig, 'Tag', 'YP1EditText' );   
   nux1EditHandle = findobj( thisFig, 'Tag', 'nux1EditText' );
   nuy1EditHandle = findobj( thisFig, 'Tag', 'nuy1EditText' ); 
   x1 = str2double( get( x1EditHandle, 'String' ) );
   y1 = str2double( get( y1EditHandle, 'String' ) );
   xp1 = str2double( get( xp1EditHandle, 'String' ) );
   yp1 = str2double( get( yp1EditHandle, 'String' ) ); 
   nux1 = str2double( get( nux1EditHandle, 'String' ) ); if isnan(nux1) nux1=0; end
   nuy1 = str2double( get( nuy1EditHandle, 'String' ) ); if isnan(nuy1) nuy1=0; end
   if( isnan(x1) || isnan(y1) || isnan(xp1) || isnan(yp1) || x1<=0 || y1<=0 )
      errordlg( 'Target condition input error!', 'ERROR', 'modal' );
      return;   
   end
   % Weights
   xwEditHandle = findobj( thisFig, 'Tag', 'XWEditText' );
   ywEditHandle = findobj( thisFig, 'Tag', 'YWEditText' );
   xpwEditHandle = findobj( thisFig, 'Tag', 'XPWEditText' );
   ypwEditHandle = findobj( thisFig, 'Tag', 'YPWEditText' );
   DwEditHandle = findobj( thisFig, 'Tag', 'DwEditText' );
   DpwEditHandle = findobj( thisFig, 'Tag', 'DpwEditText' );  
   nuxwEditHandle = findobj( thisFig, 'Tag', 'nuxwEditText' );
   nuywEditHandle = findobj( thisFig, 'Tag', 'nuywEditText' );  
   betawEditHandle = findobj( thisFig, 'Tag', 'betawEditText' );
   refwEditHandle = findobj( thisFig, 'Tag', 'refwEditText' ); 
   xw = str2double( get( xwEditHandle, 'String' ) );
   yw = str2double( get( ywEditHandle, 'String' ) );
   xpw = str2double( get( xpwEditHandle, 'String' ) );
   ypw = str2double( get( ypwEditHandle, 'String' ) );
   Dw = str2double( get( DwEditHandle, 'String' ) ); if isnan(Dw) Dw=0; end
   Dpw = str2double( get( DpwEditHandle, 'String' ) ); if isnan(Dpw) Dpw=0; end
   nuxw = str2double( get( nuxwEditHandle, 'String' ) );if isnan(nuxw) nuxw=0; end
   nuyw = str2double( get( nuywEditHandle, 'String' ) );  if isnan(nuyw) nuyw=0; end
   betaw = str2double( get( betawEditHandle, 'String' ) );if isnan(betaw) betaw=0; end
   refw = str2double( get( refwEditHandle, 'String' ) ); if isnan(refw) refw=0; end
   if( isnan(xw) || isnan(yw) || isnan(xpw) || isnan(ypw) || xw<0 || yw<0 || xpw<0 || ypw<0 )
      errordlg( 'Weights input error', 'ERROR', 'modal' );
      return;   
   end
   % Iterations
   maxIterEditHandle = findobj( thisFig, 'Tag', 'MaxIterEditText' );
   tolFunEditHandle = findobj( thisFig, 'Tag', 'TolFunEditText' );
   maxIter = str2double( get( maxIterEditHandle, 'String' ) );
   tolFun = str2double( get( tolFunEditHandle, 'String' ) );
   if( isnan(maxIter) || isnan(tolFun) || maxIter<=0 || tolFun<=0 )
      errordlg( 'Iterations input error', 'ERROR', 'modal' );
      return;   
   end
   % Find the shared userdata
%   mainFigHandle = get( thisFig, 'UserData' );
%   usrdata = get( mainFigHandle, 'UserData' );
   % -- targets
   guim.usrdata.target.x1 = x1;
   guim.usrdata.target.y1 = y1;
   guim.usrdata.target.xp1 = xp1;
   guim.usrdata.target.yp1 = yp1;
   % -- optional targets
   guim.usrdata.target.nux1 = nux1;
   guim.usrdata.target.nuy1 = nuy1;
   % -- weights
   guim.usrdata.weights.xw = xw;
   guim.usrdata.weights.yw = yw;
   guim.usrdata.weights.xpw = xpw;
   guim.usrdata.weights.ypw = ypw;
   % -- optional weights
   guim.usrdata.weights.Dw = Dw; 
   guim.usrdata.weights.Dpw = Dpw; 
   guim.usrdata.weights.nuxw = nuxw; 
   guim.usrdata.weights.nuyw = nuyw; 
   guim.usrdata.weights.betaw = betaw;
   guim.usrdata.weights.refw = refw;
   % -- matcher settings
   guim.usrdata.maxIter = maxIter;
   guim.usrdata.tolFun = tolFun;   
%   set( mainFigHandle, 'UserData', usrdata );
case 'Cancel'
end
close( thisFig );


function updateListboxString( listboxHandle, usrdata, focus )
% Add formated string to the caller listbox
[~,n] = size( usrdata.loc );
if( n==0 )
   listboxString = [];
else
   listboxString = cell(n,1);
end
for i=1:n
   if( usrdata.ele(i)=='Q' ) ele='Quad';
   elseif( usrdata.ele(i)=='S' ) ele='Sol';
   else ele='Dipl'; end
   if( usrdata.opt(i) ) ele=[ele '*']; end
   itemString = sprintf('%-7s%-8.3f  %-8.3f  %-12.6f %-8.3f',ele,usrdata.loc(i),usrdata.len(i),usrdata.str(i),usrdata.irho(i) );
   listboxString(i) = { itemString };
end
set( listboxHandle, 'String', listboxString );
set( listboxHandle, 'Value', focus );


function changeButtonState( fig, nListboxItems )
edtButtonHandle = findobj( fig, 'Tag', 'EditPushbutton' );
delButtonHandle = findobj( fig, 'Tag', 'DeletePushbutton' );
if( nListboxItems )
   set( edtButtonHandle, 'Enable', 'on' );
   set( delButtonHandle, 'Enable', 'on' );
else
   set( edtButtonHandle, 'Enable', 'off' );
   set( delButtonHandle, 'Enable', 'off' );
end

function axesHandle = findAxes( fig )
axesHandle = findobj( fig, 'Type', 'axes' );
for i=1:length(axesHandle)
   pos = get( axesHandle(i), 'Position' );
   x(i) = pos(1);
end
[x,ind] = sort(x);
% Sort axes handle according to their x position
axesHandle = axesHandle( ind );

function clearAxes( axesHandle )
axdata = get( axesHandle, 'UserData' );
for i=1:4
   if( axdata.handle(i)~=0 )
      delete( axdata.handle(i) );
      axdata.handle(i) = 0;
   end
end
set( axesHandle, 'UserData', axdata );

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
   'irho',[], ...
   'opt',[], ...
   'x1',[], ...
   'y1',[], ...
   'xp1',[], ...
   'yp1',[], ...   
   'xw',1.0, ...
   'yw',1.0, ...
   'xpw',1.0, ...
   'ypw',1.0, ...
   'Dw',0.0, ...
   'Dpw',0.0, ...
   'nuxw',0.0, ...
   'nuyw',0.0, ...
   'betaw',0.0, ...
   'refw',0.0, ...
   'xref',[], ...
   'yref',[], ...
   'maxIter',20, ...
   'tolFun',1.e-8 );

function usrdata = zeroAxesUserData
usrdata = struct( ...
   'handle',[0,0,0,0,0,0], ...
   'd',[], ...
   'x',[], ...
   'y',[] );

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
    if isfield(data.target,'D1')
        usrdata.target.D1 = data.target.D1*1.e2;           % m->cm
    end
end
% -- reference trajectory
if isfield(data.target,'xref')
    usrdata.target.xref = data.target.xref*1.e2;               % m->cm
    usrdata.target.yref = data.target.yref*1.e2;               % m->cm
end
%usrdata.xp1,usrdata.yp1,usrdata.xw,usrdata.yw,usrdata.xpw,usrdata.ypw,usrdata.maxIter,usrdata.tolFun


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
try
    usrdata.d0 = data.d0*1.e-2;   % cm->m
catch
    usrdata.d0 = 0;
end
%usrdata.ele
usrdata.loc = data.loc*1.e-2;             % cm->m
usrdata.len = data.len*1.e-2;             % cm->m
%usrdata.str
if isfield(data,'target')
    usrdata.target.x1 = data.target.x1*1.e-2;               % cm->m
    usrdata.target.y1 = data.target.y1*1.e-2;               % cm->m
    if isfield(data.target,'D1')
        usrdata.target.D1 = data.target.D1*1.e-2;           % m->cm
    end
end
if isfield(data.target,'xref')
    usrdata.target.xref = data.target.xref*1.e-2;               % m->cm
    usrdata.target.yref = data.target.yref*1.e-2;               % m->cm
end
%usrdata.xp1,usrdata.yp1,usrdata.xw,usrdata.yw,usrdata.xpw,usrdata.ypw,usrdata.maxIter,usrdata.tolFun


function savefile( usrdata, filename )
save 'savetmp' usrdata;
fcopyfile( 'savetmp.mat', filename );
delete( 'savetmp.mat' );

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
