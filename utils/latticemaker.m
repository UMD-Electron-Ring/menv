% Build Lattice; some default builds specified, included parser for MAD
% files

classdef latticemaker
    methods
        
        function [lat] = parsemad(self,madfile)
        	% -- grab info from .lat file   
            a = importdata('BTF_lattice.lat');
            accel = struct();
            for line=1:length(a)
                thisline = strsplit(a{line},{' ',',',':',';','(',')'});
                % -- first two entries are name and type
                elename = thisline{1};
                eletype = thisline{2};
                params = struct('type',eletype);
                % -- additional entries are element parameters (if any)
                if length(thisline)>2 && ~strncmpi(eletype,'LINE',4)
                   for iparam = 3:length(thisline)
                       if ~isempty(thisline{iparam}) % ignore an empty cell
                           thisparam = strsplit(thisline{iparam},{'='});
                           paramname = thisparam{1};
                           paramval = str2num(thisparam{2});
                           params = setfield(params,paramname,paramval);
                       end
                   end
                end
                if strncmpi(eletype,'LINE',4)
                    beamline = strsplit(a{line},{'=','(',')',','});
                    bl=beamline(2:end-1);
                    bl=reshape(bl,length(bl),1);
                    params = setfield(params,'beamline',bl);
                end
                accel = setfield(accel,elename,params);
            end
            lat = accel;
            
            % -- Convert to menv arrays
            
            %
            % INCOMPLETE
            %
            
            
        end
        
        
        function [lat,distance] = parsexml(self,infile,beamlinename)
            
            dat= textread(infile,'%s','endofline','\n');
            
            % -- ask for beamline name if not specified
            if not(exist('beamlinename'))
                dat % display imported file contents
                beamlinename = input('Type name of beamline to be imported from XML file: ','s');
            end
                
            
            % -- isolate just beamlinename
            istart = find(~cellfun(@isempty,strfind(dat,sprintf('<%s',beamlinename))));
            iend = find(~cellfun(@isempty,strfind(dat,sprintf('%s>',beamlinename))));
            dat = dat(istart:iend);
            
            % -- grab beamline length
            istart = find(~cellfun(@isempty,strfind(dat,sprintf('<%s',beamlinename))));
            iend = find(~cellfun(@isempty,strfind(dat,'>')),1);
            subdat = dat(istart:iend);
            ii = find(~cellfun(@isempty,strfind(subdat,'length')),1);
            bookends = strfind(subdat{ii},'"');
            beamlinelength = str2num(subdat{ii}(bookends(1)+1:bookends(2)-1));
            
            
            % -- find start/end indices for individual element definitions
            istart = find(~cellfun(@isempty,strfind(dat,'<accElement')));
            iend = find(~cellfun(@isempty,strfind(dat,'accElement>')));
            
            % -- make structure to store element values
            Nele = length(istart);
            accel = struct();
            
            for i=1:Nele
                % -- make temporary structure to save attributes
                tempstruct = struct();
                % -- select data for element i
                subdat = dat(istart(i):iend(i));
                % -- name
                ii = find(~cellfun(@isempty,strfind(subdat,'name')),1);
                bookends = strfind(subdat{ii},'"');
                name = subdat{ii}(bookends(1)+1:bookends(2)-1);
                tempstruct=setfield(tempstruct,'name',name);
                % -- Type
                ii = find(~cellfun(@isempty,strfind(subdat,'type')),1);
                bookends = strfind(subdat{ii},'"');
                type = subdat{ii}(bookends(1)+1:bookends(2)-1);
                tempstruct=setfield(tempstruct,'type',type);
                % -- length
                ii = find(~cellfun(@isempty,strfind(subdat,'length')),1);
                bookends = strfind(subdat{ii},'"');
                len = str2num(subdat{ii}(bookends(1)+1:bookends(2)-1));
                tempstruct=setfield(tempstruct,'len',len);
                % -- position
                ii = find(~cellfun(@isempty,strfind(subdat,'pos')),1);
                bookends = strfind(subdat{ii},'"');
                pos = str2num(subdat{ii}(bookends(1)+1:bookends(2)-1));
                tempstruct=setfield(tempstruct,'pos',pos);
                % -- other parameters (optional)
                ii = find(~cellfun(@isempty,strfind(subdat,'parameters')),1);
                jj = find(~cellfun(@isempty,strfind(subdat,'accElement>')),1);
                subsubdat = subdat(ii+1:jj-1);
                for nparam = 1:length(subsubdat)
                    bookends = strfind(subsubdat{nparam},'"');
                    varname = subsubdat{nparam}(1:bookends(1)-2);
                    varvalue = str2num(subsubdat{nparam}(bookends(1)+1:bookends(2)-1));
                    tempstruct=setfield(tempstruct,varname,varvalue);
                end
                
                % -- assign values to accelerator structure
                elementName = sprintf('ele%i',i);
                accel = setfield(accel,elementName,tempstruct);
            end
            
            
            % -- convert accel structure into lat arrays for menv
            ele ='';
            loc = [];
            len = [];
            str =[];
            irho = [];
            counter = 1;
            
            for i=1:Nele
                elementName = sprintf('ele%i',i);
                thisele = getfield(accel,elementName);
                % -- type
                eletype = thisele.type;
                if strncmpi(eletype,'QUAD',4) elename = 'Q';
                elseif strncmpi(eletype,'BEND',4) elename = 'D';
                else warning(sprintf('Element type %s unexpected. Skipping it, if important add to latticemaker.parsexml',eletype))
                    continue
                end
                ele=[ele,elename];
                % -- loc
                eleloc = thisele.pos * 1e2; % convert to cm
                loc(counter) = eleloc;
                % -- len
                elelen = thisele.len * 1e2; % convert to cm
                len(counter) = elelen;
                % -- str
                try
                    elestr = thisele.field; % geometric strength?
                catch elestr = 0.;
                    warning(sprintf('Focusing not included for element %s',thisele.name))
                end
                str(counter) = elestr;
                % -- irho
                try
                    eleang = thisele.theta; % bending angle
                catch eleang = 0.;
                end
                irho(counter) = eleang/elelen;
                
                counter = counter+1; % count the added element
            end
            % -- sort list
            [loc,isort] = sort(loc);
            
            len = len(isort);
            str = str(isort);
            opt = 0*len;
            irho = irho(isort);
            
            lat.ele = ele;
            lat.len = len;
            lat.loc = loc;
            lat.str = str;
            lat.opt = opt;
            lat.irho = irho;
            
            distance = beamlinelength;
        end
                
                
        
        
        function [lat] = build(ncells)
           
            Iquad = 1.826;
            nD = ncells;
            nQ = 2*ncells;
            nele = nD+nQ;
            lcell = 32;
            L = lcell*ncells;
            
            % -- dipo params
            dlen = 3.819;
            dang = 10*pi/180;
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
                opt(i) = 0;
            end
            
            % -- define dipoles
            for i=1:nD
                loc(i+nQ) = 16 + 32*(i-1);
                len(i+nQ) = dlen;
                str(i+nQ) = dstr;
                opt(i+nQ) = 0;
                irho(i+nQ) = 1/rho;
            end
            
            % -- sort list
            [loc,isort] = sort(loc);
            
            len = len(isort);
            str = str(isort);
            opt = opt(isort);
            irho = irho(isort);
            
            lat.ele = ele;
            lat.len = len;
            lat.loc = loc;
            lat.str = str;
            lat.opt = opt;
            lat.irho = irho;
            
        end
        
    end
end