function [X]=MOmatching()

% -- Load otimization settings
load 'runtmp';
maxIter = round(runtmp.maxIter);
tolFun = runtmp.tolFun;

% -- make lattmp variable and load initial opt. parameters into X0
X0 = menv_makekappaarray();

% make initial population
NPOP = 20;
% pop0 = [59.59	-51.493	-42.11	52.735];
% IP = ((rand(NPOP,length(pop0))*2)-1)*300; IP(1,:) = pop0;


options = gaoptimset('CreationFcn', @gacreationlinearfeasible,...
                    'MutationFcn', @mutationadaptfeasible,...
                    'CrossoverFcn', @crossoverintermediate,...
                    'PopulationSize',NPOP,...
                    'ParetoFraction',0.7,...
                    'PlotFcns',@gaplotpareto,...
                    'PopInitRange',[-300;300]);
%                    'InitialPopulation',pop0);
 NVARS = length(OPT_ELE);
 UB = 300*ones(1,NVARS); LB=-UB;
[x fval flag output population] = gamultiobj(@stepfuncMO,NVARS,...
                          [],[],[],[],LB,UB,options);


runtmp.f = fval;
save 'runtmp' runtmp;
                      
end


