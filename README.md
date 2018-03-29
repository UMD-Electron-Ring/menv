# menv

Matching code written by Hui Li. Based off SPOT by Chris Allen 
(does not make use of reference trajectory in optimization function)

~ User Guide ~
To run, fig = menv launches GUI. Can open .spt file, or build lattice 
inside GUI using add/edit elements.

menv.m -- Initiates GUI
menv.mat -- holds data needed to set up GUI figures.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DEV-NOTES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Original MENV was a bit disorganized -- identical routines were repeated in
multiple places, sub-function hierarchy was not consistent. Many desired
features were not implemented (dispersion, chromaticity, matching to reference trajectory)

Branch nlu-matching has been dedicated to updating MENV, starting with a command-line clone clmenv.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3/29/18

command-line menv (clmenv) is almost fully re-organized ad de-bugged. 
File hierarchy is below:
clmenv
--match2period
--match2target
--GSmatch2target
--MOmatch2target
---runmenv
---periodfunc
---optfunc
---svoptfunc
----menv_makekappaarray
----menv_integrator
----menv_updatekappaarray
----Kappa2Current
----Current2Kappa
----Gauss2Kappa
----DrawReferenceTrajectory

example scripts: run_clmenv and run_clmenv_inj

GUI files:
menv.mat
menv.m
menvEvent.m
prepareMenu.m
defElement.m
defMatcher.m
defParam.m



To-do:

--benchmark w/ old trace3D model.

--Multi-objective (MO) matching method is not yet operational. 
Needs to be finished in 'new organization'
Soon-to-be deprecated files:
MOmatching.m,stepfuncMO.m

--Implement new organization and features in MENV GUI.

--Once GUI MENV is updated I'll be comfortable deleting some more of the old, deprecated files:
step.m
calc_prim2.m
stepfunc.matching

--Incorporate skew terms: this was apparently started by a precessor but never finished?
See skew_test.m, skewstep.m. The physics included here would need to be incorporated 
into the updated organization