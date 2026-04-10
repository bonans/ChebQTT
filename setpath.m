function setpath()
%SETPATH Setup paths and check for required toolboxes.
%   SETPATH() adds the source files and external toolboxes to the MATLAB
%   path, and verifies that all required dependencies are available.
%
%   IMPORTANT: Modify the paths below to match your installation locations.
%
%   See also: REPRODUCE_ALL, RUN_COMPARE_PLOT, RUN_COMPARE_TABLE, RUN_3D_PLOTS

%% Path configuration
% =========================================================================
% IMPORTANT: Modify these paths to match your installation locations
% =========================================================================
path_to_chebfun = '~/Documents/MATLAB/chebfun/';
path_to_TT_Toolbox = '~/Documents/MATLAB/TT-Toolbox/';

%% Add source files
% Get the directory where setpath.m is located
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'src'));

%% Add external toolboxes
addpath(path_to_chebfun);
addpath(path_to_TT_Toolbox);
addpath([path_to_TT_Toolbox 'core/']);
addpath([path_to_TT_Toolbox 'exp/']);
addpath([path_to_TT_Toolbox 'cross/']);
addpath([path_to_TT_Toolbox 'fmex/']);
addpath([path_to_TT_Toolbox 'misc/']);
addpath([path_to_TT_Toolbox 'solve/']);

%% Check for required toolboxes
fprintf('Checking for required toolboxes...\n');

% Check Chebfun
if ~exist('chebfun', 'file')
    error(['Chebfun not found. Please install it from https://www.chebfun.org/ ' ...
           'and update the path_to_chebfun variable in setpath.m.']);
end
fprintf('  Chebfun: OK\n');

% Check TT-Toolbox
if ~exist('tt_tensor', 'file')
    error(['TT-Toolbox not found. Please install it from ' ...
           'https://github.com/oseledets/TT-Toolbox ' ...
           'and update the path_to_TT_Toolbox variable in setpath.m.']);
end
fprintf('  TT-Toolbox: OK\n');

%% Compile MEX files (if needed)
src_dir = fullfile(this_dir, 'src');
mex_sources = {'clenshaw_eval_mex.c', 'horner_eval_mex.c'};

all_compiled = true;
for i = 1:numel(mex_sources)
    [~, name] = fileparts(mex_sources{i});
    if ~exist(fullfile(src_dir, [name '.' mexext]), 'file')
        all_compiled = false;
        break;
    end
end

if all_compiled
    fprintf('  MEX files: OK (already compiled)\n');
else
    fprintf('  Compiling MEX files...\n');
    for i = 1:numel(mex_sources)
        src_file = fullfile(src_dir, mex_sources{i});
        fprintf('    %s\n', mex_sources{i});
        try
            mex('COPTIMFLAGS=-O3 -march=native -ffast-math -DNDEBUG', ...
                '-outdir', src_dir, src_file);
        catch ME
            error(['Failed to compile ' mex_sources{i} ':\n' ME.message]);
        end
    end
    fprintf('  MEX files: OK\n');
end

fprintf('\n');
end
