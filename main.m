%% ODE ,odel, constant cytoplasmic Plk4 concentration

clear
load('Cycle12_ODE_Oscillator_ConstantPlk4.mat')

Plk4cycle(xPlk4, t, 'multCycle', 'off', 'offset', 'off', 'plot','on', 'Plk4shape', Plk4shape, 'Plk4stage', Plk4stage, 'PPTshape', PPTshape, 'PPTstage', PPTstage, 'Sfrac', Sfrac, 'tperiod', tperiod);

%% ODE model, variable cytoplasmic Plk4 concentration

clear
load('Cycle12_ODE_Oscillator_VaryPlk4.mat')

Plk4cycle(xPlk4, t, 'multCycle', 'on', 'offset', 'off', 'plot','on', 'Plk4shape', Plk4shape, 'Plk4stage', Plk4stage, 'PPTshape', PPTshape, 'PPTstage', PPTstage, 'Sfrac', Sfrac, 'tperiod', tperiod);

%% Stochastic model, non-cycling (may take some time)

clear
load('Cycle12_Stochastic_FreeRun.mat')

Plk4cycleStochastic(xPlk4, t, 'multCycle', 'off', 'plot','on', 'tperiod', tperiod);
