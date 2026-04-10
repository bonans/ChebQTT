function run_compare_table(varargin)
%RUN_COMPARE_TABLE Generate comparison tables for joint QTT approximation.
%   RUN_COMPARE_TABLE() generates comparison tables for joint QTT
%   approximation of multiple polynomials.
%
%   Optional name-value pairs:
%     'mode', M           - Execution mode (default: 'compute&plot')
%                           'compute'      : Only compute, save data to data/
%                           'plot'         : Only print tables/CSV from saved
%                                            data (falls back to 'compute&plot'
%                                            if no data)
%                           'compute&plot' : Compute, save data, and print
%     'repeat', N         - Number of timing repeats (default: 1)
%     'saveresults', bool - Save CSV to disk (default: false)
%                           WARNING: Setting to true will overwrite existing files!
%     'd', D              - QTT dimension, grid size 2^d (default: 20)
%
%   This function generates Tables 2-3 of the paper, comparing:
%     - Separate vs Joint (First) vs Joint (Last) configurations
%     - Constructive low-rank method vs dmrg_cross
%
%   For test functions:
%     1. Synthetic: 50 random Chebyshev polynomials of degree 200
%     2. Biomolecule potential (ChebTuck format)
%     3. Runge function (ChebTuck format)
%     4. Wagon function (ChebTuck format)
%
%   See also: RUN_COMPARE_PLOT, RUN_3D_PLOTS, SETPATH

% Parse inputs
p = inputParser;
addParameter(p, 'mode', 'compute&plot', @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
addParameter(p, 'repeat', 1, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'saveresults', false, @islogical);
addParameter(p, 'd', 20, @(x) isnumeric(x) && x >= 1);
parse(p, varargin{:});
mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});
repeat = p.Results.repeat;
saveresults = p.Results.saveresults;
d = p.Results.d;

% Determine what to do
datafile = fullfile('data', 'compare_table_results.mat');
do_compute = contains(mode, 'compute');
do_plot = contains(mode, 'plot');

% Fallback: if plot-only but no data, switch to compute&plot
if strcmp(mode, 'plot') && ~isfile(datafile)
    warning('No saved data found (%s). Falling back to compute&plot mode.', datafile);
    do_compute = true;
end

% In plot-only mode, load parameters from saved data
if ~do_compute && do_plot && isfile(datafile)
    saved = load(datafile, 'd');
    d = saved.d;
end

if do_compute
    fprintf('=== Running joint QTT comparison tables ===\n');
else
    fprintf('=== Skipped computation of joint QTT comparison tables (plotting from saved data) ===\n');
end
fprintf('Parameters: mode=%s, repeat=%d, saveresults=%s, d=%d\n\n', ...
    mode, repeat, mat2str(saveresults), d);

t_wall_start = tic;

%% Parameters
methods = {'ConstructiveLowRank', 'TTCross'};
method_labels = {'Constructive', 'dmrg_cross'};
n_methods = length(methods);
configs = {'Separate', 'JointFirst', 'JointLast'};
config_labels = {'Separate', 'Joint First', 'Joint Last'};
n_configs = length(configs);
n_functions = 4;     % Synthetic, Biomolecule, Runge, Wagon
func_names = {'Synthetic', 'Biomolecule', 'Runge', 'Wagon'};

%% Computation
if do_compute

% Results: 24 rows x 4 columns (error, storage, avg rank, time)
% Row indexing: (func_idx-1)*6 + (config_idx-1)*2 + method_idx
results = zeros(n_functions * n_configs * n_methods, 4);

%% Function loop
for func_idx = 1:n_functions
    fprintf('=== Function %d/%d: %s ===\n', func_idx, n_functions, func_names{func_idx});

    % Load/generate polynomial coefficients
    switch func_idx
        case 1  % Synthetic: 50 random Chebyshev polynomials of degree 200
            rng(114514);
            M = 50; m = 200;
            a = 0.1*rand(1,M); b = 10*rand(1,M);
            funs_mean = b.*exp(-a.*(0:m)');
            funs_coef = funs_mean + 0.1.*funs_mean.*randn(m+1,M);

        case 2  % Biomolecule potential
            % Load native data (Chebyshev coefficients) and reconstruct chebfun3
            biomol_data = load("data/biomol.mat");
            f = chebfun3();
            f.domain = biomol_data.domain;
            f.cols = chebfun(biomol_data.cols_coeffs, biomol_data.domain(1:2), 'coeffs');
            f.rows = chebfun(biomol_data.rows_coeffs, biomol_data.domain(3:4), 'coeffs');
            f.tubes = chebfun(biomol_data.tubes_coeffs, biomol_data.domain(5:6), 'coeffs');
            f.core = biomol_data.core;
            funs_coef = biomol_data.cols_coeffs;

        case 3  % Runge function
            ff = @(x,y,z) 1./(1 + 25*(x.^2 + y.^2 + z.^2));
            f = chebfun3(ff);
            funs_coef = chebcoeffs(f.cols);

        case 4  % Wagon function
            ff = @(x,y,z) exp(sin(50*x)) + sin(60*exp(y)).*sin(60*z) + ...
                sin(70*sin(x)).*cos(10*z) + sin(sin(80*y)) - sin(10*(x+z)) ...
                + (x.^2 + y.^2 + z.^2)/4;
            f = chebfun3(ff);
            funs_coef = chebcoeffs(f.cols);
    end

    M = size(funs_coef, 2);
    funs_coef_cell = {funs_coef, 'cheb'};
    fprintf('  %d polynomials of degree %d\n', M, size(funs_coef,1)-1);

    dom = [-1, 1];

    % Print table header
    fprintf('%-12s %-15s %15s %10s %10s %12s\n', ...
        'Method', 'Configuration', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
    fprintf('%s\n', repmat('-', 1, 76));

    %% Configuration and method loop (config outer → methods batched)
    for config_idx = 1:n_configs
        config = configs{config_idx};

        err_vals     = zeros(n_methods, 1);
        storage_vals = zeros(n_methods, 1);
        rank_vals    = zeros(n_methods, 1);
        time_vals    = zeros(n_methods, 1);
        method_ok    = true(n_methods, 1);

        switch config
            case 'Separate'
                % Build all M TTs for each method
                all_tts = cell(n_methods, M);
                for method_idx = 1:n_methods
                    method = methods{method_idx};
                    for jj = 1:M
                        try
                            tt = poly2qtt({funs_coef(:,jj), 'cheb'}, d, 'method', method, 'order', 'first');
                            all_tts{method_idx, jj} = tt;
                            storage_vals(method_idx) = storage_vals(method_idx) + mem(tt);
                            rank_vals(method_idx) = rank_vals(method_idx) + TTrankMean(tt)/M;
                        catch ME
                            fprintf('      Warning: Failed for poly %d: %s\n', jj, ME.message);
                            method_ok(method_idx) = false;
                            err_vals(method_idx) = NaN;
                            storage_vals(method_idx) = NaN;
                            rank_vals(method_idx) = NaN;
                            break;
                        end
                    end
                    % Timing
                    for r = 1:repeat
                        t_start = tic;
                        for jj = 1:M
                            try
                                poly2qtt({funs_coef(:,jj), 'cheb'}, d, 'method', method, 'order', 'first');
                            catch
                            end
                        end
                        time_vals(method_idx) = time_vals(method_idx) + toc(t_start)/repeat;
                    end
                end

                % Batch error per column (ground truth evaluated once per column)
                for jj = 1:M
                    col_ok_idx = find(method_ok & ~cellfun(@isempty, all_tts(:, jj)));
                    if ~isempty(col_ok_idx)
                        errs = qtt_poly_error(all_tts(col_ok_idx, jj), funs_coef(:,jj), 'cheb', dom);
                        for k = 1:numel(col_ok_idx)
                            err_vals(col_ok_idx(k)) = max(err_vals(col_ok_idx(k)), errs(k));
                        end
                    end
                end

            case {'JointFirst', 'JointLast'}
                if strcmp(config, 'JointFirst')
                    order_str = 'first';
                else
                    order_str = 'last';
                end

                tts = cell(n_methods, 1);
                for method_idx = 1:n_methods
                    method = methods{method_idx};
                    try
                        tt = poly2qtt(funs_coef_cell, d, 'method', method, 'order', order_str);
                        tts{method_idx} = tt;
                        storage_vals(method_idx) = mem(tt);
                        rank_vals(method_idx) = TTrankMean(tt);
                    catch ME
                        fprintf('      Warning: Failed: %s\n', ME.message);
                        method_ok(method_idx) = false;
                        err_vals(method_idx) = NaN;
                        storage_vals(method_idx) = NaN;
                        rank_vals(method_idx) = NaN;
                    end
                    % Timing
                    for r = 1:repeat
                        t_start = tic;
                        try
                            poly2qtt(funs_coef_cell, d, 'method', method, 'order', order_str);
                        catch
                        end
                        time_vals(method_idx) = time_vals(method_idx) + toc(t_start)/repeat;
                    end
                end

                % Batch error computation (ground truth evaluated once)
                ok_idx = find(method_ok);
                if ~isempty(ok_idx)
                    errs = qtt_poly_error(tts(ok_idx), funs_coef, 'cheb', dom, 'order', order_str);
                    for k = 1:numel(ok_idx)
                        err_vals(ok_idx(k)) = errs(k);
                    end
                end
        end

        % Store results and print
        for method_idx = 1:n_methods
            row_idx = (func_idx-1)*6 + (config_idx-1)*2 + method_idx;
            results(row_idx, :) = [err_vals(method_idx), storage_vals(method_idx), ...
                rank_vals(method_idx), time_vals(method_idx)];
            fprintf('%-12s %-15s %15.2e %10d %10.2f %12.4f\n', ...
                method_labels{method_idx}, config_labels{config_idx}, ...
                err_vals(method_idx), storage_vals(method_idx), ...
                rank_vals(method_idx), time_vals(method_idx));
        end
    end
end

% Save computed data
if ~exist('data', 'dir'), mkdir('data'); end
save(datafile, 'results', 'methods', 'method_labels', 'configs', 'config_labels', ...
    'func_names', 'n_methods', 'n_configs', 'n_functions', 'd', 'repeat');
fprintf('Data saved to: %s\n', datafile);

end % do_compute

%% Print formatted tables
if do_plot

% Load data if not already in workspace (plot-only mode)
if ~do_compute
    fprintf('Loading data from: %s\n', datafile);
    loaded = load(datafile);
    results = loaded.results;
    methods = loaded.methods;
    method_labels = loaded.method_labels;
    config_labels = loaded.config_labels;
    func_names = loaded.func_names;
    n_methods = loaded.n_methods;
    n_configs = loaded.n_configs;
    n_functions = loaded.n_functions;
    d = loaded.d;
end

fprintf('\n');
fprintf('=========================================================================\n');
fprintf('                    JOINT QTT APPROXIMATION RESULTS                     \n');
fprintf('=========================================================================\n');

for func_idx = 1:n_functions
    fprintf('\n--- %s ---\n', func_names{func_idx});
    fprintf('%-12s %-15s %15s %10s %10s %12s\n', ...
        'Method', 'Configuration', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
    fprintf('%s\n', repmat('-', 1, 76));

    for config_idx = 1:n_configs
        for method_idx = 1:n_methods
            row_idx = (func_idx-1)*6 + (config_idx-1)*2 + method_idx;
            fprintf('%-12s %-15s %15.2e %10d %10.2f %12.4f\n', ...
                method_labels{method_idx}, config_labels{config_idx}, ...
                results(row_idx, 1), results(row_idx, 2), ...
                results(row_idx, 3), results(row_idx, 4));
        end
    end
end

%% Save results to CSV
if saveresults
    figdir = 'figures/';
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    row_names = cell(24, 1);
    method_short = {'Constr', 'dmrg_cross'};
    config_short = {'Sep', 'JFirst', 'JLast'};
    for func_idx = 1:n_functions
        for config_idx = 1:n_configs
            for method_idx = 1:n_methods
                row_idx = (func_idx-1)*6 + (config_idx-1)*2 + method_idx;
                row_names{row_idx} = sprintf('%s_%s_%s', func_names{func_idx}, config_short{config_idx}, method_short{method_idx});
            end
        end
    end

    col_names = {'Error', 'Storage', 'AvgRank', 'Time'};
    filename = fullfile(figdir, 'joint_qtt_results.csv');
    matrix2csv(results, filename, col_names, row_names);
    fprintf('\nResults saved to: %s\n', filename);
end

end % do_plot

t_wall = toc(t_wall_start);
fprintf('\n=== Joint QTT comparison complete (%.1f seconds / %.1f minutes) ===\n', t_wall, t_wall/60);
end
