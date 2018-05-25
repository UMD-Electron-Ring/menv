function menvEvent( event )
switch( get(gcbf,'Tag') )
case 'MainFig'
   MainFigWndProc( gcbf, event );
   clm = clmenv;
case 'DefElementFig'
   DefElementWndProc( gcbf, event );
case 'DefParamFig'
   DefParamWndProc( gcbf, event );
case 'DefMatcherFig'
   DefMatcherWndProc( gcbf, event );
end

function MainFigWndProc( thisFig, event )
switch( event )
case 'Create'
   usrdata = zeroMainUserData;
   set( thisFig, 'UserData', usrdata );
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
    [filename, pathname] = uigetfile('*.spt', 'Open File');
%    if( filename==0 ) return; end
%    len = length( filename );
%    if( len<5 || ~strcmpi(filename(len-3:len),'.spt') )
%       errordlg( 'The filename must have an extension .spt!', 'Error', 'modal');
%       return;
%    end
%    filename = [pathname filename];
%    fcopyfile( filename, 'savetmp.mat' );
%    load 'savetmp';
%    usrdata.filename = filename;  % important to rename its filename
%    % To be compatible with previous version file with no usrdata.did
%    if( isfield(usrdata,'did')==0 )
%       usrdata.did = zeros( size(usrdata.opt) );
%    end
%    % Transfer from SI
%    usrdata = TransferFromSI( usrdata );
   clm.open(filename,pathname)
   set( thisFig, 'UserData', clm.usrdata );
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, usrdata, 1 );
   [~,n] = size( usrdata.loc );
   changeButtonState( thisFig, n );
case 'Close'
   % Clear current axis
   axesHandle = findAxes( thisFig );
   clearAxes( axesHandle(1) );
   % Clear userData
   usrdata = zeroMainUserData;
   set( thisFig, 'UserData', usrdata );
   % Clear listbox
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   set( listboxHandle , 'String', '' );
   changeButtonState( thisFig, 0 );
case 'Save'
   usrdata = get( thisFig, 'UserData' );
   if( isempty(usrdata.filename) )
      menvEvent( 'SaveAs' );
   else
      % Transfer to SI
      usrdata = Transfer2SI( usrdata );
      savefile( usrdata, usrdata.filename );
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
   usrdata.filename = filename;
   set( thisFig, 'UserData', usrdata );
   % Transfer to SI
   usrdata = Transfer2SI( usrdata );   
   savefile( usrdata, filename );
case 'Edit'
   usrdata = get( thisFig, 'UserData' );
   % Pass the calling instruction   
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   item = get( listboxHandle, 'Value' );
   usrdata.flag = item;
   set( thisFig, 'UserData', usrdata );
   % Call the input dialogbox
   fig = defElement;
   set( fig, 'UserData', thisFig );
   % Fill out the input dialogbox
   elePopupmenuHandle = findobj( fig, 'Tag', 'ElementPopupMenu' );
   if( usrdata.ele(item)=='Q' )
      set( elePopupmenuHandle, 'Value', 1 );
   elseif( usrdata.ele(item)=='S' )
      set( elePopupmenuHandle, 'Value', 2 );
   else
      set( elePopupmenuHandle, 'Value', 3 );
   end
   locEditHandle = findobj( fig, 'Tag', 'LocationEditText' );
   lenEditHandle = findobj( fig, 'Tag', 'LengthEditText' );
   strEditHandle = findobj( fig, 'Tag', 'StrengthEditText' );
   didEditHandle = findobj( fig, 'Tag', 'DiplIndexEditText' );
   locString = num2str( usrdata.loc(item) );
   lenString = num2str( usrdata.len(item) );
   strString = num2str( usrdata.str(item) );
   didString = num2str( usrdata.did(item) );
   set( locEditHandle, 'String', locString );
   set( lenEditHandle, 'String', lenString ); 
   set( strEditHandle, 'String', strString );
   set( didEditHandle, 'String', didString );
   optCheckboxHandle = findobj( fig, 'Tag', 'OptimCheckbox' );
   set( optCheckboxHandle, 'Value', usrdata.opt(item) );
   % important in here to disable/enable some controls
   DefElementWndProc( fig, 'ElementChange' ); 
case 'Insert'
   usrdata = get( thisFig, 'UserData' );
   usrdata.flag = 0;
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
   usrdata.did = usrdata.did(ind);
   usrdata.opt = usrdata.opt(ind);
   set( thisFig, 'UserData', usrdata );
   % Update the listbox
   updateListboxString( listboxHandle, usrdata, 1 );
   changeButtonState( thisFig, n-1 );
case 'Param'
   fig = defParam;
   set( fig, 'UserData', thisFig );
   usrdata = get( thisFig, 'UserData' );
   % Initial conditions
   if( ~isempty(usrdata.x0) )
      x0EditHandle = findobj( fig, 'Tag', 'X0EditText' );
      y0EditHandle = findobj( fig, 'Tag', 'Y0EditText' );
      xp0EditHandle = findobj( fig, 'Tag', 'XP0EditText' );
      yp0EditHandle = findobj( fig, 'Tag', 'YP0EditText' );
      set( x0EditHandle, 'String', num2str(usrdata.x0) );
      set( y0EditHandle, 'String', num2str(usrdata.y0) );
      set( xp0EditHandle, 'String', num2str(usrdata.xp0) );
      set( yp0EditHandle, 'String', num2str(usrdata.yp0) );
   end
   % Global beam parameters
   if( ~isempty(usrdata.emitance) )
      emitEditHandle = findobj( fig, 'Tag', 'EmitEditText' );
      set( emitEditHandle, 'String', num2str(usrdata.emitance) );
   end
   if( ~isempty(usrdata.perveance) )
      pvnEditHandle = findobj( fig, 'Tag', 'PvnEditText' );
      set( pvnEditHandle, 'String', num2str(usrdata.perveance) );
   end
   % Numerical parameters
   if( ~isempty(usrdata.stepsize) )
      stepEditHandle = findobj( fig, 'Tag', 'StepEditText' );
      set( stepEditHandle, 'String', num2str(usrdata.stepsize) );
   end
   if( ~isempty(usrdata.distance) )
      distEditHandle = findobj( fig, 'Tag', 'DistEditText' );
      set( distEditHandle, 'String', num2str(usrdata.distance) );
   end
case 'Solution'
   runtmp = get( thisFig, 'UserData' );
   % only check one parameter is enough
   if( isempty(runtmp.x0) )
      warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   % Transfer to SI
   runtmp = Transfer2SI( runtmp );
   % Save to tempary file
   save 'runtmp' runtmp;
   % Run ......
   oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
   [x,y,xp,yp,D,Dp,d,nux,nuy,Cx,Cy] = runmenv( 'runtmp' );
   x=x*1.e2; y=y*1.e2; d=d*1.e2;  % m-cm
   D=D*1.e2;
   set( thisFig, oldptr{:} );
   % Axes
   axesHandle = findAxes( thisFig );
   axes( axesHandle(1) );
   axdata = get( axesHandle(1), 'UserData' );
   % Plot
   if( axdata.handle(1)~=0 ) delete(axdata.handle(1)); end
   if( axdata.handle(2)~=0 ) delete(axdata.handle(2)); end
   hold on; h1 = plot(d,x,'b'); h2 = plot(d,y,'r'); hold off;
   xlabel('z (cm)'); ylabel('X:blue, Y:red (cm)');
   axis([ min(d) max(d) 0.0 max([x,y])*1.2 ]);
   
    % Save data
   if( usrdata.stepsize<0.005 ) interval = round(0.005/usrdata.stepsize);
   else interval = 1; end
   n = length(d); ind = 1:interval:n;
   if( ind(length(ind))~=n ) ind(length(ind)+1) = n; end   
   axdata.handle(1)=h1; axdata.handle(2)=h2; 
   axdata.d=d(ind); axdata.x=x(ind); axdata.y=y(ind);
   axdata.xp=xp(ind); axdata.yp=yp(ind);
   axdata.nux = nux; axdata.nuy = nuy;
   set( axesHandle, 'UserData', axdata );   
case 'MatcherParam'
   fig = defMatcher;
   set( fig, 'UserData', thisFig );
   usrdata = get( thisFig, 'UserData' );
   % Target conditions
   if( ~isempty(usrdata.x1) )
      x1EditHandle = findobj( fig, 'Tag', 'X1EditText' );
      y1EditHandle = findobj( fig, 'Tag', 'Y1EditText' );
      xp1EditHandle = findobj( fig, 'Tag', 'XP1EditText' );
      yp1EditHandle = findobj( fig, 'Tag', 'YP1EditText' );
      set( x1EditHandle, 'String', num2str(usrdata.x1) );
      set( y1EditHandle, 'String', num2str(usrdata.y1) );
      set( xp1EditHandle, 'String', num2str(usrdata.xp1) );
      set( yp1EditHandle, 'String', num2str(usrdata.yp1) );
   end
   % Weights
   if( ~isempty(usrdata.xw) )
      xwEditHandle = findobj( fig, 'Tag', 'XWEditText' );
      ywEditHandle = findobj( fig, 'Tag', 'YWEditText' );
      xpwEditHandle = findobj( fig, 'Tag', 'XPWEditText' );
      ypwEditHandle = findobj( fig, 'Tag', 'YPWEditText' );
      set( xwEditHandle, 'String', num2str(usrdata.xw) );
      set( ywEditHandle, 'String', num2str(usrdata.yw) );
      set( xpwEditHandle, 'String', num2str(usrdata.xpw) );
      set( ypwEditHandle, 'String', num2str(usrdata.ypw) );
   end
   % iterations
   if( ~isempty(usrdata.maxIter) )
      maxIterEditHandle = findobj( fig, 'Tag', 'MaxIterEditText' );
      tolFunEditHandle = findobj( fig, 'Tag', 'TolFunEditText' );
      set( maxIterEditHandle, 'String', num2str(usrdata.maxIter) );
      set( tolFunEditHandle, 'String', num2str(usrdata.tolFun) );
   end
case 'Run Periodic Matcher'
   usrdata = get( thisFig, 'UserData' );
   % only check one parameter is enough
   if( isempty(usrdata.x0) )
      warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   if( isempty(usrdata.x1) )
      warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   % Transfer to SI
   usrdata = Transfer2SI( usrdata );
   % Save to tempary file
   save 'runtmp' usrdata;
   % Run ......
   %oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
   oldptr = get(thisFig,'pointer');  set( thisFig, 'pointer', 'watch' );
   newX0 = match2period( 'runtmp' );
   set( thisFig, 'pointer', oldptr );
   % Save the new result
   usrdata.x0 = newX0(1);
   usrdata.y0 = newX0(2);
   usrdata.xp0= newX0(3);
   usrdata.yp0= newX0(4);
   usrdata = TransferFromSI( usrdata );
   set( thisFig, 'UserData', usrdata );
   % Update figure
   menvEvent( 'Solution' );
case 'Matcher'
   usrdata = get( thisFig, 'UserData' );
   % only check one parameter is enough
   if( isempty(usrdata.x0) )
      warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   if( isempty(usrdata.x1) )
      warndlg( 'Matcher Parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   if( sum(usrdata.opt)==0 )
      warndlg( 'The optimized elements are not defined!', 'ERROR', 'modal' );
      return;
   end
   % Transfer to SI
   usrdata = Transfer2SI( usrdata );
   % Save to tempary file
   save 'runtmp' usrdata;
   % Run ......
   oldptr = getptr(thisFig);  setptr( thisFig, 'watch' );
   newKappa = match2target( 'runtmp' );
   set( thisFig, oldptr{:} );
   % Save the new result
   usrdata = get( thisFig, 'UserData' );
   [~,n] = size( usrdata.loc ); k = 1;
   for i=1:n
      if( usrdata.opt(i) )
         usrdata.str(i) = newKappa(k); k = k+1;
      end
   end
   set( thisFig, 'UserData', usrdata );
   % Update the listbox
   listboxHandle = findobj( thisFig, 'Tag', 'ElementListbox' );
   updateListboxString( listboxHandle, usrdata, 1 );   
   % Update figure
   menvEvent( 'Solution' );
case 'Coordinate'
   coordTextHandle = findobj( thisFig, 'Tag', 'CoordStaticText' );
   axesHandle = findAxes( thisFig );
   xlim = get( axesHandle, 'XLim' );
   ylim = get( axesHandle, 'YLim' );
   set( coordTextHandle, 'Visible', 'on' );
   while( 1 )
      [x,y] = ginput(1);
      str = sprintf('(%f, %f)',x,y);
      set( coordTextHandle, 'String', str );
      if( x<xlim(1) || x>xlim(2) || y<ylim(1) || y>ylim(2) ) break; end
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
   usrdata = get( thisFig, 'UserData' );
   % only check one parameter is enough
   if( isempty(usrdata.x0) )
      warndlg( 'Beam parameters are not defined!', 'ERROR', 'modal' );
      return;
   end
   % Transfer to SI
   usrdata = Transfer2SI( usrdata );
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
    usrdata = get( thisFig, 'UserData' );
    % Draw on figure
    axesHandle = findAxes( thisFig );
    [xref,yref] = DrawReferenceTrajectory(axesHandle);
    % save in UserData
    usrdata.xref = xref;
    usrdata.yref = yref;
    set( thisFig, 'UserData', usrdata );
    % Save new reference trajectory to file
    [filename, pathname] = uiputfile('*.txt', 'Save As');
    dlmwrite([pathname,filename],[xref,yref]);
    % Plot on axes
    hold on
    if exist('href'), delete(href); end
    href = plot(xref,yref,':k');
    hold off
    
case 'Load Reference Trajectory'
    usrdata = get( thisFig, 'UserData' );    
    % load .txt file
    [filename, pathname] = uigetfile('*.txt', 'Open File');
    if( filename==0 ) return; end
    len = length( filename );
    filename = [pathname filename];
    temp =  importdata(filename);
    xref = temp(:,1); yref = temp(:,2);    
    % save in UserData
    usrdata.xref = xref;
    usrdata.yref = yref;
    set( thisFig, 'UserData', usrdata );
    % Plot on axes
    hold on
    if exist('href'), delete(href); end
    href = plot(xref,yref,':k');
    hold off
    
otherwise
end


function DefElementWndProc( thisFig, event )
switch( event )
case 'OK'
   % Judge if all the inputs valid
   elePopupmenuHandle = findobj( thisFig, 'Tag', 'ElementPopupMenu' );
   locEditHandle = findobj( thisFig, 'Tag', 'LocationEditText' );
   lenEditHandle = findobj( thisFig, 'Tag', 'LengthEditText' );
   strEditHandle = findobj( thisFig, 'Tag', 'StrengthEditText' );
   didEditHandle = findobj( thisFig, 'Tag', 'DiplIndexEditText' );
   element = get( elePopupmenuHandle, 'Value' );
   locString = get( locEditHandle, 'String' );
   lenString = get( lenEditHandle, 'String' );
   strString = get( strEditHandle, 'String' );
   didString = get( didEditHandle, 'String' );
   location = str2double( locString );
   length = str2double( lenString );
   strength = str2double( strString );
   diplindex = str2double( didString );
   if( isnan(location) || isnan(length) || isnan(strength) || (element==3 && isnan(diplindex)) )
      errordlg( 'Some input parameters are not numbers!', 'ERROR', 'modal' );
      return;
   end
   optCheckboxHandle = findobj( thisFig, 'Tag', 'OptimCheckbox' );
   optim = get( optCheckboxHandle, 'Value' );
   % Find the shared userdata
   mainFigHandle = get( thisFig, 'UserData' );
   listboxHandle = findobj( mainFigHandle, 'Tag', 'ElementListbox' );
   usrdata = get( mainFigHandle, 'UserData' );
   % 1:Quad; 2:Sol; 3:Dipl
   if( element==1 )
      element = 'Q';
      diplindex = 0;
   elseif( element==2 )
      element = 'S';
      diplindex = 0;
   else
      element = 'D';
      optim = 0; % disable optim for Dipl
   end
   % 0:Insert; 1:Edit
   if( usrdata.flag==0 )
      usrdata.ele = [usrdata.ele element];
      usrdata.loc = [usrdata.loc location];
      usrdata.len = [usrdata.len length];
      usrdata.str = [usrdata.str strength];
      usrdata.did = [usrdata.did diplindex];
      usrdata.opt = [usrdata.opt optim];
   else
      item = usrdata.flag;
      usrdata.ele(item) = element;
      usrdata.loc(item) = location;
      usrdata.len(item) = length;
      usrdata.str(item) = strength;
      usrdata.did(item) = diplindex;
      usrdata.opt(item) = optim;
   end
   % Find who is the item that we are editing or inserting
   [~,n] = size( usrdata.loc ); 
   [tmp,ind] = sort( usrdata.loc );
   if( usrdata.flag==0 )
      focus = find( ind==n );
   else
      focus = find( ind==item );
   end
   % Sort
   usrdata.loc = tmp; 
   usrdata.ele = usrdata.ele(ind);
   usrdata.len = usrdata.len(ind);
   usrdata.str = usrdata.str(ind);
   usrdata.did = usrdata.did(ind);
   usrdata.opt = usrdata.opt(ind);
   % Save
   set( mainFigHandle, 'UserData', usrdata );
   % Add formated string to the caller listbox
   updateListboxString( listboxHandle, usrdata, focus );
   changeButtonState( mainFigHandle, n );
   close( thisFig );
case 'Cancel'
   close( thisFig );
case 'ElementChange'
   elePopupmenuHandle = findobj( thisFig, 'Tag', 'ElementPopupMenu' );
   didEditHandle = findobj( thisFig, 'Tag', 'DiplIndexEditText' );
   optCheckboxHandle = findobj( thisFig, 'Tag', 'OptimCheckbox' );
   element = get( elePopupmenuHandle, 'Value' );
   if( element==3 ) % Dipl
      set( optCheckboxHandle, 'Enable', 'off' );
      set( didEditHandle, 'Enable', 'on' );
      set( didEditHandle, 'BackgroundColor', [1 1 1] ); % white
   else
      set( optCheckboxHandle, 'Enable', 'on' );
      set( didEditHandle, 'Enable', 'off' );
      set( didEditHandle, 'BackgroundColor', [0.7529 0.7529 0.7529] ); % gray      
   end
end


function DefParamWndProc( thisFig, event )
switch( event )
case 'OK'
   % Initial conditions
   x0EditHandle = findobj( thisFig, 'Tag', 'X0EditText' );
   y0EditHandle = findobj( thisFig, 'Tag', 'Y0EditText' );
   xp0EditHandle = findobj( thisFig, 'Tag', 'XP0EditText' );
   yp0EditHandle = findobj( thisFig, 'Tag', 'YP0EditText' );   
   x0 = str2double( get( x0EditHandle, 'String' ) );
   y0 = str2double( get( y0EditHandle, 'String' ) );
   xp0 = str2double( get( xp0EditHandle, 'String' ) );
   yp0 = str2double( get( yp0EditHandle, 'String' ) );   
   if( isnan(x0) || isnan(y0) || isnan(xp0) || isnan(yp0) || x0<=0 || y0<=0 )
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
   mainFigHandle = get( thisFig, 'UserData' );
   usrdata = get( mainFigHandle, 'UserData' );
   usrdata.emitance = emitance;
   usrdata.perveance = perveance;
   usrdata.x0 = x0;
   usrdata.y0 = y0;
   usrdata.xp0 = xp0;
   usrdata.yp0 = yp0;
   usrdata.stepsize = stepsize;
   usrdata.distance = distance;
   set( mainFigHandle, 'UserData', usrdata );
case 'Cancel'
end
close( thisFig );

function DefMatcherWndProc( thisFig, event )
switch( event )
case 'OK'
   % Target conditions
   x1EditHandle = findobj( thisFig, 'Tag', 'X1EditText' );
   y1EditHandle = findobj( thisFig, 'Tag', 'Y1EditText' );
   xp1EditHandle = findobj( thisFig, 'Tag', 'XP1EditText' );
   yp1EditHandle = findobj( thisFig, 'Tag', 'YP1EditText' );   
   x1 = str2double( get( x1EditHandle, 'String' ) );
   y1 = str2double( get( y1EditHandle, 'String' ) );
   xp1 = str2double( get( xp1EditHandle, 'String' ) );
   yp1 = str2double( get( yp1EditHandle, 'String' ) );   
   if( isnan(x1) || isnan(y1) || isnan(xp1) || isnan(yp1) || x1<=0 || y1<=0 )
      errordlg( 'Target condition input error!', 'ERROR', 'modal' );
      return;   
   end
   % Weights
   xwEditHandle = findobj( thisFig, 'Tag', 'XWEditText' );
   ywEditHandle = findobj( thisFig, 'Tag', 'YWEditText' );
   xpwEditHandle = findobj( thisFig, 'Tag', 'XPWEditText' );
   ypwEditHandle = findobj( thisFig, 'Tag', 'YPWEditText' );   
   xw = str2double( get( xwEditHandle, 'String' ) );
   yw = str2double( get( ywEditHandle, 'String' ) );
   xpw = str2double( get( xpwEditHandle, 'String' ) );
   ypw = str2double( get( ypwEditHandle, 'String' ) );   
   if( isnan(xw) || isnan(yw) || isnan(xpw) || isnan(ypw) || xw<=0 || yw<=0 || xpw<=0 || ypw<=0 )
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
   mainFigHandle = get( thisFig, 'UserData' );
   usrdata = get( mainFigHandle, 'UserData' );
   usrdata.x1 = x1;
   usrdata.y1 = y1;
   usrdata.xp1 = xp1;
   usrdata.yp1 = yp1;
   usrdata.xw = xw;
   usrdata.yw = yw;
   usrdata.xpw = xpw;
   usrdata.ypw = ypw;
   usrdata.maxIter = maxIter;
   usrdata.tolFun = tolFun;   
   set( mainFigHandle, 'UserData', usrdata );
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
   itemString = sprintf('%-7s%-8.3f  %-8.3f  %-12.6f',ele,usrdata.loc(i),usrdata.len(i),usrdata.str(i) );
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
   'tolFun',1.e-8 );

function usrdata = zeroAxesUserData
usrdata = struct( ...
   'handle',[0,0,0,0,0,0], ...
   'd',[], ...
   'x',[], ...
   'y',[] );

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
