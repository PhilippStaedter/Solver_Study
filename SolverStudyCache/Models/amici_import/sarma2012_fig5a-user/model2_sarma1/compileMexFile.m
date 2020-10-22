% This compile script was automatically created from Python SBML import.
% If mex compiler is set up within MATLAB, it can be run from MATLAB 
% in order to compile a mex-file from the Python generated C++ files.

modelName = 'model2_sarma1';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(27, 27, 45, 0, 0, 0, [], ...
            ['simulate_' modelName '.m'], modelName, 'lin', 1, 1);