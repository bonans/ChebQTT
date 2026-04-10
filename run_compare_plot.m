function run_compare_plot(varargin)
%RUN_COMPARE_PLOT Compare QTT approximation methods for single polynomials.
%   RUN_COMPARE_PLOT() generates comparison plots for different QTT
%   approximation methods across various polynomial types.
%
%   Optional name-value pairs:
%     'mode', M           - Execution mode (default: 'compute&plot')
%                           'compute'      : Only compute, save data to data/
%                           'plot'         : Only plot from saved data (falls
%                                            back to 'compute&plot' if no data)
%                           'compute&plot' : Compute, save data, and plot
%     'repeat', N         - Number of timing repeats (default: 1)
%     'saveresults', bool - Save figures to disk (default: false)
%                           WARNING: Setting to true will overwrite existing files!
%     'd', D              - QTT dimension, grid size 2^d (default: 20)
%     'm_max', M          - Maximum polynomial degree (default: 300)
%
%   This function generates Figure 3 of the paper, showing approximation
%   errors, average QTT ranks, and runtimes for:
%     (A) Monomials
%     (B) Random polynomials in monomial basis
%     (C) Chebyshev polynomials
%     (D) Random polynomials in Chebyshev basis
%
%   See also: RUN_COMPARE_TABLE, RUN_3D_PLOTS, SETPATH

% Parse inputs
p = inputParser;
addParameter(p, 'mode', 'compute&plot', @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
addParameter(p, 'repeat', 1, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'saveresults', false, @islogical);
addParameter(p, 'd', 20, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'm_max', 300, @(x) isnumeric(x) && x >= 1);
parse(p, varargin{:});
mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});
repeat = p.Results.repeat;
saveresults = p.Results.saveresults;
d = p.Results.d;
m_max = p.Results.m_max;

% Determine what to do
datafile = fullfile('data', 'compare_plot_results.mat');
do_compute = contains(mode, 'compute');
do_plot = contains(mode, 'plot');

% Fallback: if plot-only but no data, switch to compute&plot
if strcmp(mode, 'plot') && ~isfile(datafile)
    warning('No saved data found (%s). Falling back to compute&plot mode.', datafile);
    do_compute = true;
end

% In plot-only mode, load parameters from saved data
if ~do_compute && do_plot && isfile(datafile)
    saved = load(datafile, 'd', 'm_max');
    d = saved.d;
    m_max = saved.m_max;
end

if do_compute
    fprintf('=== Running single polynomial comparison ===\n');
else
    fprintf('=== Skipped computation of single polynomial comparison (plotting from saved data) ===\n');
end
fprintf('Parameters: mode=%s, repeat=%d, saveresults=%s, d=%d, m_max=%d\n\n', ...
    mode, repeat, mat2str(saveresults), d, m_max);

t_wall_start = tic;

methods = {'ConstructiveMonomial', 'ConstructiveLowRank', ...
    'ConstructiveFullRank', 'TTCross', 'Horner'};
method_labels = {'Constructive (mono)', 'Constructive (low-rank)', ...
    'Constructive (full-rank)', 'TT-Cross', 'Horner/Clenshaw'};
n_methods = length(methods);
% Testing sets of polynomials
% single mono, random mono, single cheb, random cheb
n_set = 4;
set_names = {'Monomials', 'Random (monomial basis)', ...
    'Chebyshev polynomials', 'Random (Chebyshev basis)'};
set_short = {'mono', 'mono_rand', 'cheb', 'cheb_rand'};

%% Computation
if do_compute

errors = zeros(m_max, n_methods, n_set);
runtimes = zeros(m_max, n_methods, n_set);
meanranks = zeros(m_max, n_methods, n_set);

% Generate random coefficients
rng(114514)
a = 0.1 * rand();
b = 10 * rand();
funs_mean = b * exp(-a * (0:m_max)');
rand_coef = funs_mean + 0.1 .* funs_mean .* randn(m_max + 1, 1);

%% Main computation loop
method_short = {'C(mono)', 'C(LR)', 'C(FR)', 'TTCross', 'Horner'};
dom = [-1, 1];

for s = 1:n_set
    is_cheb_set = (s >= 3);
    n_active = n_methods - is_cheb_set;

    % Print table header for this set
    fprintf('\n--- Set %d/%d: %s (d=%d) ---\n', s, n_set, set_names{s}, d);
    fprintf('%6s', 'm');
    for i = 1:n_methods
        if is_cheb_set && strcmpi(methods{i}, 'ConstructiveMonomial'), continue; end
        fprintf('  %12s', method_short{i});
    end
    fprintf('\n%s\n', repmat('-', 1, 6 + n_active * 14));

    for m = 1:m_max
        % Setup polynomial coefficients
        if s == 1 || s == 3  % single monomial or Chebyshev
            funs_coef = zeros(m + 1, 1);
            funs_coef(end) = 1;
            funs_coef = {funs_coef, 'mono'};
            if s == 3
                funs_coef = {funs_coef{1}, 'cheb'};
            end
        else  % random polynomials
            funs_coef = {rand_coef(1:m + 1), 'mono'};
            if s == 4  % random Chebyshev
                funs_coef = {rand_coef(1:m + 1), 'cheb'};
            end
        end

        % Collect TTs for all methods, then batch-compute errors
        tts = cell(n_methods, 1);
        method_ok = false(n_methods, 1);
        for i = 1:n_methods
            method = methods{i};

            % Skip Chebyshev + ConstructiveMonomial (not applicable)
            if strcmpi(funs_coef{2}, 'cheb') && strcmpi(method, 'ConstructiveMonomial')
                continue;
            end

            % First run (compile/warmup)
            tt_method = poly2qtt(funs_coef, d, 'method', method);
            tts{i} = tt_method;
            method_ok(i) = true;
            meanranks(m, i, s) = TTrankMean(tt_method);

            % Timing runs
            for r = 1:repeat
                tic
                tt_method = poly2qtt(funs_coef, d, 'method', method);
                runtimes(m, i, s) = runtimes(m, i, s) + toc;
            end
        end

        % Batch error computation (ground truth evaluated once)
        ok_idx = find(method_ok);
        if ~isempty(ok_idx)
            errs = qtt_poly_error(tts(ok_idx), funs_coef{1}, funs_coef{2}, dom);
            for k = 1:numel(ok_idx)
                errors(m, ok_idx(k), s) = errs(k);
            end
        end

        % Print progress row (L∞ errors for each method)
        fprintf('%6d', m);
        for i = 1:n_methods
            if is_cheb_set && strcmpi(methods{i}, 'ConstructiveMonomial'), continue; end
            if errors(m, i, s) > 0 || meanranks(m, i, s) > 0
                fprintf('  %12.2e', errors(m, i, s));
            else
                fprintf('  %12s', '---');
            end
        end
        fprintf('\n');
    end
end

% Average the runtimes
runtimes = runtimes / repeat;

% Save computed data
if ~exist('data', 'dir'), mkdir('data'); end
save(datafile, 'errors', 'runtimes', 'meanranks', ...
    'methods', 'method_labels', 'set_names', 'set_short', ...
    'd', 'm_max', 'repeat', 'n_methods', 'n_set');
fprintf('Data saved to: %s\n', datafile);

end % do_compute

%% Plotting
if do_plot

% Load data if not already in workspace (plot-only mode)
if ~do_compute
    fprintf('Loading data from: %s\n', datafile);
    loaded = load(datafile);
    errors = loaded.errors;
    runtimes = loaded.runtimes;
    meanranks = loaded.meanranks;
    methods = loaded.methods;
    method_labels = loaded.method_labels;
    set_names = loaded.set_names;
    set_short = loaded.set_short;
    d = loaded.d;
    m_max = loaded.m_max;
    n_methods = loaded.n_methods;
    n_set = loaded.n_set;
end

% Julia color scheme (matching paper figures)
julia_colors = [
    0.0, 0.6056031611752245, 0.978680117569607;      % julia1: blue (TTCross)
    0.8888735002725198, 0.43564919034818983, 0.2781229361419438;  % julia2: orange (ConstructiveMonomial)
    0.2422242978521988, 0.6432750931576304, 0.30444865153411527;  % julia3: green (ConstructiveFullRank)
    0.7644401754934356, 0.44411177946877645, 0.824297535923276;   % julia4: purple (ConstructiveLowRank)
    0.6755439572114058, 0.5556623322045814, 0.09423433626639477;  % julia5: yellow-brown (Horner)
];

% Method order in 'methods' array:
% 1: ConstructiveMonomial -> julia2 (orange)
% 2: ConstructiveLowRank  -> julia4 (purple)
% 3: ConstructiveFullRank -> julia3 (green)
% 4: TTCross              -> julia1 (blue)
% 5: Horner               -> julia5 (yellow-brown)
method_colors = [julia_colors(2,:); julia_colors(4,:); julia_colors(3,:); julia_colors(1,:); julia_colors(5,:)];

% Y-axis limits for each set (matching paper TikZ singlepoly.tex)
% Set 1: mono, Set 2: mono_rand, Set 3: cheb, Set 4: cheb_rand
error_ylim = {[1e-14, 1e-4], [1e-14, 1e-4], [1e-14, 1e-4], [1e-14, 1e-7]};
rank_ylim = {[1, 8], [1, 8], [1, 17], [1, 17]};
runtime_ylim = {[0, 0.5], [0, 0.5], [0, 0.5], [0, 0.5]};

figdir = 'figures/';
if saveresults && ~exist(figdir, 'dir')
    mkdir(figdir);
end

% Combined figure: 4 columns (A-D) x 3 rows (error, rank, runtime)
% Layout matches paper TikZ: top section = mono (A,B), bottom = cheb (C,D)
% We use a 6x2 subplot grid: rows 1-3 = mono (A left, B right),
%                              rows 4-6 = cheb (C left, D right)
fig = figure('Position', [100, 100, 800, 900]);

% Column order: [mono, mono_rand, cheb, cheb_rand] -> subplot columns [1,2,1,2]
% Row offset: mono sets (s=1,2) -> rows 1-3, cheb sets (s=3,4) -> rows 4-6
col_map = [1, 2, 1, 2];       % subplot column for each set
row_off = [0, 0, 3, 3];       % row offset for each set
n_rows = 6; n_cols = 2;

for s = 1:n_set
    is_cheb = (s >= 3);
    col = col_map(s);
    roff = row_off(s);

    % Row 1: Error
    subplot(n_rows, n_cols, (roff)*n_cols + col);
    hold on;
    for i = 1:n_methods
        if is_cheb && strcmpi(methods{i}, 'ConstructiveMonomial'), continue; end
        valid_idx = errors(:, i, s) > 0;
        if any(valid_idx)
            semilogy(find(valid_idx), errors(valid_idx, i, s), ...
                'Color', method_colors(i,:), 'LineWidth', 1.5, 'DisplayName', method_labels{i});
        end
    end
    hold off;
    set(gca, 'YScale', 'log', 'FontSize', 9);
    ylim(error_ylim{s});
    if s == 4
        set(gca, 'YTick', 10.^(-14:2:-7));
    else
        set(gca, 'YTick', 10.^(-14:2:-4));
    end
    xlim([1, m_max]);
    set(gca, 'XTickLabel', []);
    grid on; box on;
    if col == 1
        ylabel('$\ell_\infty$ error', 'Interpreter', 'latex', 'FontSize', 10);
    end
    % Title only for row 1
    if roff == 0
        title(sprintf('\\textbf{(%s)} %s', char('A' + s - 1), set_names{s}), ...
            'Interpreter', 'latex', 'FontSize', 10);
    else
        title(sprintf('\\textbf{(%s)} %s', char('A' + s - 1), set_names{s}), ...
            'Interpreter', 'latex', 'FontSize', 10);
    end
    if s == 1
        legend('Location', 'best', 'FontSize', 7);
    end

    % Row 2: Rank
    subplot(n_rows, n_cols, (roff+1)*n_cols + col);
    hold on;
    for i = 1:n_methods
        if is_cheb && strcmpi(methods{i}, 'ConstructiveMonomial'), continue; end
        valid_idx = meanranks(:, i, s) > 0;
        if any(valid_idx)
            plot(find(valid_idx), meanranks(valid_idx, i, s), ...
                'Color', method_colors(i,:), 'LineWidth', 1.5);
        end
    end
    hold off;
    set(gca, 'FontSize', 9);
    ylim(rank_ylim{s});
    if s <= 2
        set(gca, 'YTick', 0:2:8);
    else
        set(gca, 'YTick', 0:4:20);
    end
    xlim([1, m_max]);
    set(gca, 'XTickLabel', []);
    grid on; box on;
    if col == 1
        ylabel('avg. rank $\bar{r}$', 'Interpreter', 'latex', 'FontSize', 10);
    end

    % Row 3: Runtime
    subplot(n_rows, n_cols, (roff+2)*n_cols + col);
    hold on;
    for i = 1:n_methods
        if is_cheb && strcmpi(methods{i}, 'ConstructiveMonomial'), continue; end
        valid_idx = runtimes(:, i, s) > 0;
        if any(valid_idx)
            plot(find(valid_idx), runtimes(valid_idx, i, s), ...
                'Color', method_colors(i,:), 'LineWidth', 1.5);
        end
    end
    hold off;
    set(gca, 'FontSize', 9);
    ylim(runtime_ylim{s});
    set(gca, 'YTick', 0:0.1:0.4);
    xlim([1, m_max]);
    grid on; box on;
    if col == 1
        ylabel('runtime (s)', 'Interpreter', 'latex', 'FontSize', 10);
    end
    % Only show x-label on bottom row of each section
    if roff == 3
        xlabel('polynomial degree $m$', 'Interpreter', 'latex', 'FontSize', 10);
    else
        set(gca, 'XTickLabel', []);
    end
end

if saveresults
    filename = fullfile(figdir, 'singlepoly.pdf');
    save_figure_tight(fig, filename);
    fprintf('Saved: %s\n', filename);
end

end % do_plot

t_wall = toc(t_wall_start);
fprintf('\n=== Single polynomial comparison complete (%.1f seconds / %.1f minutes) ===\n', t_wall, t_wall/60);
end

function save_figure_tight(fig_handle, filename)
    fig_handle.PaperUnits = 'inches';
    fig_pos = fig_handle.Position;
    screen_dpi = 100;
    fig_width_inches = fig_pos(3) / screen_dpi;
    fig_height_inches = fig_pos(4) / screen_dpi;
    fig_handle.PaperSize = [fig_width_inches, fig_height_inches];
    fig_handle.PaperPosition = [0, 0, fig_width_inches, fig_height_inches];
    print(fig_handle, filename, '-dpdf', '-bestfit');
end
