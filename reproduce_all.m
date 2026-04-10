function reproduce_all(varargin)
%REPRODUCE_ALL Reproduce all figures and tables from the paper.
%   REPRODUCE_ALL() runs all experiments with default parameters for fast
%   prototyping (~4 minutes).
%
%   REPRODUCE_ALL('paper') uses the exact parameters from the paper
%   (d=20, m_max=300, repeat=20, saveresults=true). Expected runtime: ~20 hours.
%
%   Optional name-value pairs:
%     'mode', M           - 'compute&plot' (default), 'compute', or 'plot'
%     'repeat', N         - Number of timing repeats (default: 1)
%     'saveresults', bool - Save figures/CSV to disk (default: false)
%     'd', D              - QTT dimension, grid size 2^d (default: 10)
%     'm_max', M          - Max polynomial degree for Experiment 1 (default: 100)
%     'd_values', [...]   - Grid dimensions for Experiment 4 (default: 10:2:20)
%
%   Presets:
%     reproduce_all()         - Fast prototyping (d=10, m_max=100, repeat=1, d_values=10:2:20)
%     reproduce_all('paper')  - Paper reproduction (d=20, m_max=300, repeat=20, d_values=10:2:30)
%
%   REQUIREMENTS:
%     - MATLAB R2022a or later
%     - Chebfun package (https://www.chebfun.org/)
%     - TT-Toolbox (https://github.com/oseledets/TT-Toolbox)
%
%   See also: RUN_COMPARE_PLOT, RUN_COMPARE_TABLE, RUN_3D_PLOTS, RUN_SCALING, SETPATH

% Check for preset shortcut: reproduce_all('paper')
use_paper = false;
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) && strcmpi(varargin{1}, 'paper')
    use_paper = true;
    varargin(1) = [];
end

% Default values (fast prototyping)
default_mode = 'compute&plot';
default_repeat = 1;
default_save = false;
default_d = 10;
default_m_max = 100;
default_d_values = 10:2:20;

% Paper values override defaults
if use_paper
    default_repeat = 20;
    default_save = true;
    default_d = 20;
    default_m_max = 300;
    default_d_values = 10:2:30;
end

% Parse inputs
p = inputParser;
addParameter(p, 'mode', default_mode, @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
addParameter(p, 'repeat', default_repeat, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'saveresults', default_save, @islogical);
addParameter(p, 'd', default_d, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'm_max', default_m_max, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'd_values', default_d_values, @isnumeric);
parse(p, varargin{:});

mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});
repeat = p.Results.repeat;
saveresults = p.Results.saveresults;
d = p.Results.d;
m_max = p.Results.m_max;
d_values = p.Results.d_values;

%% Setup paths
setpath();

% In plot-only mode, load parameters from saved data for correct display
if strcmp(mode, 'plot')
    cp_file = fullfile('data', 'compare_plot_results.mat');
    if isfile(cp_file)
        saved = load(cp_file, 'd', 'm_max');
        d = saved.d;
        m_max = saved.m_max;
    end
end

%% Run experiments
fprintf('=========================================================================\n');
fprintf('   Reproducing results for:\n');
fprintf('   "Bridging continuous and discrete tensor formats of multivariate\n');
fprintf('    functions via QTT"\n');
fprintf('=========================================================================\n\n');
fprintf('Configuration:\n');
fprintf('  mode       = %s\n', mode);
fprintf('  repeat     = %d\n', repeat);
fprintf('  saveresults = %s\n', mat2str(saveresults));
fprintf('  d          = %d (grid size 2^d = %d)\n', d, 2^d);
fprintf('  m_max      = %d\n', m_max);
fprintf('  d_values   = %s (Experiment 4)\n', mat2str(d_values));
if strcmp(mode, 'plot')
    fprintf('  (parameters loaded from saved data)\n');
end
fprintf('\n');

total_time = tic;

% Snapshot pre-existing files in figures/ for summary
figdir = 'figures/';
pre_existing_files = {};
if exist(figdir, 'dir')
    pre_files = [dir(fullfile(figdir, '*.pdf')); dir(fullfile(figdir, '*.csv'))];
    pre_existing_files = {pre_files.name};
end

%% Experiment 1: Single polynomial comparison (Figure 3)
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('   Experiment 1: Single polynomial QTT approximation (Figure 3)\n');
fprintf('=========================================================================\n');
fprintf('This experiment compares different methods for constructing QTT\n');
fprintf('representations of single polynomials up to degree %d.\n\n', m_max);

run_compare_plot('mode', mode, 'repeat', repeat, 'saveresults', saveresults, 'd', d, 'm_max', m_max);
if saveresults
    close all;
end

%% Experiment 2: Joint QTT comparison tables (Tables 2-3)
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('   Experiment 2: Joint QTT approximation (Tables 2-3)\n');
fprintf('=========================================================================\n');
fprintf('This experiment compares separate vs joint QTT approximation\n');
fprintf('for polynomial sets from Chebyshev-Tucker decompositions.\n\n');

run_compare_table('mode', mode, 'repeat', repeat, 'saveresults', saveresults, 'd', d);
if saveresults
    close all;
end

%% Experiment 3: 3D surface visualizations (Figure 4)
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('   Experiment 3: 3D surface visualization (Figure 4)\n');
fprintf('=========================================================================\n');
fprintf('This experiment visualizes the QTT-Tucker approximations of\n');
fprintf('multivariate functions on a %dx%d grid.\n\n', 2^min(d-1,10), 2^min(d-1,10));

run_3d_plots('mode', mode, 'saveresults', saveresults, 'd', d);
if saveresults
    close all;
end

%% Experiment 4: Scaling with respect to d (Figure 5)
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('   Experiment 4: Scaling with respect to d (Figure 5)\n');
fprintf('=========================================================================\n');
fprintf('This experiment tests scaling behavior as the quantized dimension d\n');
fprintf('varies over d = %s for random Chebyshev polynomials.\n\n', mat2str(d_values));

run_scaling('mode', mode, 'repeat', repeat, 'saveresults', saveresults, 'd_values', d_values);
if saveresults
    close all;
end

%% Summary
elapsed = toc(total_time);
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('   COMPLETE\n');
fprintf('=========================================================================\n');
fprintf('Total runtime: %.1f seconds (%.1f minutes)\n\n', elapsed, elapsed/60);

% Report file status in figures/
figdir = 'figures/';
if saveresults && exist(figdir, 'dir')
    post_files = [dir(fullfile(figdir, '*.pdf')); dir(fullfile(figdir, '*.csv'))];
    if isempty(post_files)
        fprintf('No files were saved to %s.\n', figdir);
    else
        fprintf('Files in %s:\n', figdir);
        for i = 1:length(post_files)
            fname = post_files(i).name;
            if ismember(fname, pre_existing_files)
                fprintf('  [replaced]   %s\n', fname);
            else
                fprintf('  [new]        %s\n', fname);
            end
        end
        % Check for pre-existing files that are still there untouched
        for i = 1:length(pre_existing_files)
            if ~ismember(pre_existing_files{i}, {post_files.name})
                fprintf('  [untouched]  %s\n', pre_existing_files{i});
            end
        end
    end
elseif ~saveresults
    if exist(figdir, 'dir')
        existing = [dir(fullfile(figdir, '*.pdf')); dir(fullfile(figdir, '*.csv'))];
        if ~isempty(existing)
            fprintf('Files in %s (not modified, saveresults=false):\n', figdir);
            for i = 1:length(existing)
                fprintf('  [untouched]  %s\n', existing(i).name);
            end
        else
            fprintf('No files in %s (nothing saved, saveresults=false).\n', figdir);
        end
    else
        fprintf('No figures/ directory (nothing saved, saveresults=false).\n');
    end
end

fprintf('\nNote: Figures in the paper were generated using TikZ/pgfplots.\n');
fprintf('This reproducibility package generates equivalent MATLAB figures.\n');
end
