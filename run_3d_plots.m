function run_3d_plots(varargin)
%RUN_3D_PLOTS Generate 3D surface plots of QTT-Tucker approximations.
%   RUN_3D_PLOTS() generates 3D visualizations of QTT-Tucker approximations
%   for multivariate functions.
%
%   Optional name-value pairs:
%     'mode', M           - Execution mode (default: 'compute&plot')
%                           'compute'      : Only compute, save data to data/
%                           'plot'         : Only plot from saved data (falls
%                                            back to 'compute&plot' if no data)
%                           'compute&plot' : Compute, save data, and plot
%     'saveresults', bool - Save figures to disk (default: false)
%                           WARNING: Setting to true will overwrite existing files!
%     'd', D              - QTT dimension, grid size 2^d (default: 20)
%
%   This function generates Figure 4 of the paper, showing:
%     - Surface plots of QTT-Tucker reconstructions (middle z-slice)
%     - Error plots comparing to continuous ChebTuck surrogates
%
%   Test functions:
%     1. Biomolecule potential (protein Fasciculin 1)
%     2. Runge function: 1/(1 + 25(x^2 + y^2 + z^2))
%     3. Wagon function: SIAM 100-Dollar, 100-Digit Challenge
%
%   See also: RUN_COMPARE_PLOT, RUN_COMPARE_TABLE, SETPATH

% Parse inputs
p = inputParser;
addParameter(p, 'mode', 'compute&plot', @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
addParameter(p, 'saveresults', false, @islogical);
addParameter(p, 'd', 20, @(x) isnumeric(x) && x >= 1);
parse(p, varargin{:});
mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});
saveresults = p.Results.saveresults;
d = p.Results.d;

% Determine what to do
datafile = fullfile('data', '3d_plots_results.mat');
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
    fprintf('=== Running 3D visualization experiments ===\n');
else
    fprintf('=== Skipped computation of 3D visualization (plotting from saved data) ===\n');
end
fprintf('Parameters: mode=%s, saveresults=%s, d=%d\n\n', mode, mat2str(saveresults), d);

t_wall_start = tic;

%% Parameters
d_plot = min(d - 1, 10);  % Plot grid dimension (smaller than d)
N_full = 2^d;             % Full grid size
N_plot = 2^d_plot;        % Plot grid size

% Subgrid: compute directly without allocating full 2^d grid
h_full = 2 / N_full;
stride = N_full / N_plot;
% Subgrid points: centered subsampling at indices stride/2, 3*stride/2, ...
uni_grid = -1 + h_full * ((stride/2 - 0.5) : stride : (N_full - stride/2 + 0.5))';
uni_grid = uni_grid(1:N_plot);  % ensure exactly N_plot points

[xgrid, ygrid] = meshgrid(uni_grid, uni_grid);
z_middle = uni_grid(N_plot/2);  % Middle slice

% Evaluation points
x_in = zeros(N_plot, N_plot, 3);
x_in(:,:,1) = xgrid;
x_in(:,:,2) = ygrid;
x_in(:,:,3) = z_middle;

func_names = {'Biomolecule', 'Runge', 'Wagon'};
func_short = {'biomol', 'runge', 'wagon'};
n_funcs = length(func_names);

%% Computation
if do_compute

% Preallocate cell arrays for surface data
surfaces = cell(n_funcs, 1);    % QTT-Tucker evaluated surfaces
references = cell(n_funcs, 1);  % Reference solutions
max_errors = zeros(n_funcs, 1);
constr_times = zeros(n_funcs, 1);
eval_times = zeros(n_funcs, 1);

% Print table header
fprintf('%-15s %12s %12s %15s\n', ...
    'Function', 'Constr (s)', 'Eval (s)', 'Max Error');
fprintf('%s\n', repmat('-', 1, 56));

%% Function 1: Biomolecule potential

% Load native data (Chebyshev coefficients) and reconstruct chebfun3
biomol_data = load("data/biomol.mat");
f1 = chebfun3();
f1.domain = biomol_data.domain;
f1.cols = chebfun(biomol_data.cols_coeffs, biomol_data.domain(1:2), 'coeffs');
f1.rows = chebfun(biomol_data.rows_coeffs, biomol_data.domain(3:4), 'coeffs');
f1.tubes = chebfun(biomol_data.tubes_coeffs, biomol_data.domain(5:6), 'coeffs');
f1.core = biomol_data.core;

tic
TuckQTT1 = ChebTuck2TuckQTT(f1, 'd', d, 'method', 'ConstructiveLowRank');
constr_times(1) = toc;

tic
surfaces{1} = TuckQTT_eval(TuckQTT1, x_in)';
eval_times(1) = toc;

references{1} = f1(xgrid, ygrid, z_middle);
max_errors(1) = max(abs(references{1}(:) - surfaces{1}(:)));
fprintf('%-15s %12.2f %12.2f %15.2e\n', func_names{1}, constr_times(1), eval_times(1), max_errors(1));

%% Function 2: Runge function
ff2 = @(x,y,z) 1./(1 + 25*(x.^2 + y.^2 + z.^2));
f2 = chebfun3(ff2);

tic
TuckQTT2 = ChebTuck2TuckQTT(f2, 'd', d, 'method', 'ConstructiveLowRank');
constr_times(2) = toc;

tic
surfaces{2} = TuckQTT_eval(TuckQTT2, x_in)';
eval_times(2) = toc;

references{2} = f2(xgrid, ygrid, z_middle);
max_errors(2) = max(abs(references{2}(:) - surfaces{2}(:)));
fprintf('%-15s %12.2f %12.2f %15.2e\n', func_names{2}, constr_times(2), eval_times(2), max_errors(2));

%% Function 3: Wagon function
ff3 = @(x,y,z) exp(sin(50*x)) + sin(60*exp(y)).*sin(60*z) + ...
    sin(70*sin(x)).*cos(10*z) + sin(sin(80*y)) - sin(10*(x+z)) ...
    + (x.^2 + y.^2 + z.^2)/4;
f3 = chebfun3(ff3);

tic
TuckQTT3 = ChebTuck2TuckQTT(f3, 'd', d, 'method', 'ConstructiveLowRank');
constr_times(3) = toc;

tic
surfaces{3} = TuckQTT_eval(TuckQTT3, x_in)';
eval_times(3) = toc;

references{3} = f3(xgrid, ygrid, z_middle);
max_errors(3) = max(abs(references{3}(:) - surfaces{3}(:)));
fprintf('%-15s %12.2f %12.2f %15.2e\n', func_names{3}, constr_times(3), eval_times(3), max_errors(3));

% Save computed data
if ~exist('data', 'dir'), mkdir('data'); end
save(datafile, 'surfaces', 'references', 'max_errors', ...
    'constr_times', 'eval_times', ...
    'uni_grid', 'func_names', 'func_short', 'n_funcs', 'd', 'd_plot');
fprintf('Data saved to: %s\n', datafile);

end % do_compute

%% Plotting
if do_plot

% Load data if not already in workspace (plot-only mode)
if ~do_compute
    fprintf('Loading data from: %s\n', datafile);
    loaded = load(datafile);
    surfaces = loaded.surfaces;
    references = loaded.references;
    max_errors = loaded.max_errors;
    if isfield(loaded, 'constr_times')
        constr_times = loaded.constr_times;
        eval_times = loaded.eval_times;
    end
    uni_grid = loaded.uni_grid;
    func_names = loaded.func_names;
    func_short = loaded.func_short;
    n_funcs = loaded.n_funcs;
    d = loaded.d;
end

% Create figures directory if needed
figdir = 'figures/';
if saveresults && ~exist(figdir, 'dir')
    mkdir(figdir);
end

for k = 1:n_funcs
    % Surface plot
    figure('Position', [100, 100, 600, 500]);
    mesh(uni_grid, uni_grid, surfaces{k})
    set(gca, 'fontsize', 16);
    view(3); grid on;
    light; lighting phong; material dull; camlight('left');
    shading interp; axis tight;
    title(sprintf('%s', func_names{k}), 'Interpreter', 'latex', 'FontSize', 18);
    if saveresults
        save_figure_tight(gcf, fullfile(figdir, [func_short{k} '.pdf']));
        fprintf('Saved: %s\n', fullfile(figdir, [func_short{k} '.pdf']));
    end

    % Error plot
    figure('Position', [100, 100, 600, 500]);
    mesh(uni_grid, uni_grid, references{k} - surfaces{k})
    set(gca, 'fontsize', 16);
    view(3); grid on;
    light; lighting phong; material dull; camlight('left');
    shading interp; axis tight;
    title(sprintf('%s error', func_names{k}), 'Interpreter', 'latex', 'FontSize', 18);
    if saveresults
        save_figure_tight(gcf, fullfile(figdir, [func_short{k} '_err.pdf']));
        fprintf('Saved: %s\n', fullfile(figdir, [func_short{k} '_err.pdf']));
    end

    fprintf('%s max error: %e\n', func_names{k}, max_errors(k));
end

% Print summary table
fprintf('\n%-15s %12s %12s %15s\n', ...
    'Function', 'Constr (s)', 'Eval (s)', 'Max Error');
fprintf('%s\n', repmat('-', 1, 56));
for k = 1:n_funcs
    if exist('constr_times', 'var') && exist('eval_times', 'var')
        fprintf('%-15s %12.2f %12.2f %15.2e\n', ...
            func_names{k}, constr_times(k), eval_times(k), max_errors(k));
    else
        fprintf('%-15s %12s %12s %15.2e\n', ...
            func_names{k}, '---', '---', max_errors(k));
    end
end

end % do_plot

t_wall = toc(t_wall_start);
fprintf('=== 3D visualization complete (%.1f seconds / %.1f minutes) ===\n', t_wall, t_wall/60);
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
