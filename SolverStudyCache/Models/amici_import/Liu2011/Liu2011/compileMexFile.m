% This compile script was automatically created from Python SBML import.
% If mex compiler is set up within MATLAB, it can be run from MATLAB 
% in order to compile a mex-file from the Python generated C++ files.

modelName = 'Liu2011';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(42, 42, 89, 0, 0, 0, [], ...
            ['simulate_' modelName '.m'], modelName, 'lin', 1, 1);