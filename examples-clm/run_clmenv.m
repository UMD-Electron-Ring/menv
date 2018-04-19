global clm
clm = clmenv();

% -- load .spt file
sptfilename = which('FODO_0_7mA_dipl.spt');
clm.open(sptfilename)

% % -- optional: makes files that may already exist
% clm.makeoptiset(20,1e-6)  % Makes file optiset
% clm.maketarget([1,1,0,0],[1,1,1,1]) % makes file target
% %clm.makeparams(); % not implemented yet, makes file params
% 
% % -- loads matcher and beam parameters into usrdata; defaults are what was
% % set in .spt file.
% clm.defparams()
% clm.defmatcher()

% -- solve
clm.periodicmatcher()

% -- save matching parameters to file
temp = clm.soldata;

%%
% -- load a different .spt file
clm.open('INJ_minsize_100mA_dipl.spt')
clm.makeoptiset(20,1e-6)
clm.maketarget([temp.xf,temp.yf,temp.xpf,temp.ypf],[1,1,1,1])
clm.defmatcher()
clm.targetmatcher()



% -- write data
nuxFODO = temp.nux;
nuyFODO = temp.nuy;
xFODO   = temp.xf;
yFODO   = temp.yf;
xpFODO  = temp.xpf;
ypFODO  = temp.ypf;
Ncells  = 10;

nuNL   = 0.23;
nuxMS   = clm.soldata.nux;
nuyMS   = clm.soldata.nuy;
nuxLIN  = nuxMS*2 + nuxFODO*Ncells - nuNL;
nuyLIN  = nuyMS*2 + nuyFODO*Ncells - nuNL;



% datatowrite = [nuxFODO, nuyFODO, Ncells*nuxFODO, Ncells*nuyFODO, ...
%     nuxMS, nuyMS, nuxLIN, nuyLIN];
% 
% fid = fopen('matching_output.csv','a');
% for i=1:length(datatowrite)
%     fprintf(fid,'%.4f,',datatowrite(i));
% end
% fprintf(fid,'\r\n');
% fclose(fid);