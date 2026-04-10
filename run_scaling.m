function run_scaling(varargin)
%RUN_SCALING Scaling experiment: QTT approximation over increasing d.
%   RUN_SCALING() sweeps d from 10 to 30 for random Chebyshev polynomials
%   of degree m=200, testing both single polynomial (M=1) and joint
%   multi-polynomial (M=50) cases. Each repeat uses a freshly generated
%   random polynomial, and results are reported as mean +/- std with error
%   bars in the plots.
%
%   Optional name-value pairs:
%     'mode', M           - 'compute', 'plot', 'compute&plot' (default)
%     'repeat', N         - Independent repeats, each with a different
%                           random polynomial (default: 20)
%     'saveresults', bool - Save figures to disk (default: false)
%     'd_values', [...]   - Grid dimensions to test (default: 10:2:30)
%     'm', M              - Polynomial degree (default: 200)
%     'M_multi', M        - Number of polynomials for joint case (default: 50)
%
%   Methods tested:
%     Single (M=1): ConstructiveLowRank, ConstructiveFullRank, TTCross, Horner
%     Joint (M=50): ConstructiveLowRank, TTCross (order='last')
%
%   Unlike run_compare_plot and run_compare_table, each repeat generates
%   a different random polynomial, so the error bars reflect variability
%   across polynomial instances.
%
%   For large d (where full(tt) would exceed available memory), the error
%   is computed via chunked evaluation (qtt_error_chunked).
%
%   See also: RUN_COMPARE_PLOT, RUN_COMPARE_TABLE, QTT_ERROR_CHUNKED

% Parse inputs
p = inputParser;
addParameter(p, 'mode', 'compute&plot', @(x) any(validatestring(x, {'compute', 'plot', 'compute&plot'})));
addParameter(p, 'repeat', 20, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'saveresults', false, @islogical);
addParameter(p, 'd_values', 10:2:30, @isnumeric);
addParameter(p, 'm', 200, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'M_multi', 50, @(x) isnumeric(x) && x >= 1);
parse(p, varargin{:});
mode = validatestring(p.Results.mode, {'compute', 'plot', 'compute&plot'});
repeat = p.Results.repeat;
saveresults = p.Results.saveresults;
d_values = p.Results.d_values;
m = p.Results.m;
M_multi = p.Results.M_multi;

% Determine what to do
datafile = fullfile('data', 'scaling_results.mat');
do_compute = contains(mode, 'compute');
do_plot = contains(mode, 'plot');

if strcmp(mode, 'plot') && ~isfile(datafile)
    warning('No saved data found (%s). Falling back to compute&plot mode.', datafile);
    do_compute = true;
end

if ~do_compute && do_plot && isfile(datafile)
    saved = load(datafile, 'd_values', 'm', 'M_multi');
    d_values = saved.d_values;
    m = saved.m;
    M_multi = saved.M_multi;
end

n_d = length(d_values);

% Methods
methods_single = {'ConstructiveLowRank', 'ConstructiveFullRank', 'TTCross', 'Horner'};
labels_single  = {'Constructive LR', 'Constructive FR', 'TT-Cross', 'Horner/Clenshaw'};
methods_joint  = {'ConstructiveLowRank', 'TTCross'};
labels_joint   = {'Constructive LR', 'TT-Cross'};
n_single = length(methods_single);
n_joint  = length(methods_joint);

fprintf('=== Scaling Experiment ===\n');
fprintf('d values: %s\n', mat2str(d_values));
fprintf('m: %d\n', m);
fprintf('M (multi): %d\n', M_multi);
fprintf('Repeat: %d, mode: %s\n\n', repeat, mode);

t_wall_start = tic;

%% Computation
if do_compute

% Results storage: dim = (d_idx, repeat_idx, method_idx)
single_error   = NaN(n_d, repeat, n_single);
single_runtime = NaN(n_d, repeat, n_single);
single_rank    = NaN(n_d, repeat, n_single);
single_storage = NaN(n_d, repeat, n_single);

joint_error   = NaN(n_d, repeat, n_joint);
joint_runtime = NaN(n_d, repeat, n_joint);
joint_rank    = NaN(n_d, repeat, n_joint);
joint_storage = NaN(n_d, repeat, n_joint);

% Fix random seed for reproducibility
rng(114514);

t_compute_start = tic;

for rr = 1:repeat
    t_repeat_start = tic;

    %% Generate one set of random polynomials per repeat (shared across all d)
    fprintf('\n========================================\n');
    fprintf('Repeat %d/%d: generating random polynomials (m=%d, M=%d)\n', ...
        rr, repeat, m, M_multi);
    fprintf('========================================\n');

    a_m = 0.1 * rand(1, M_multi);
    b_m = 10 * rand(1, M_multi);
    means_m = b_m .* exp(-a_m .* (0:m)');
    coef_multi = means_m + 0.1 .* means_m .* randn(m+1, M_multi);
    coef_single = coef_multi(:, 1);

    for di = 1:n_d
        d = d_values(di);

        %% ===== SINGLE POLYNOMIAL (M=1) =====
        funs_coef = {coef_single, 'cheb'};
        tts_single = cell(n_single, 1);
        single_ok = false(n_single, 1);
        for si = 1:n_single
            method = methods_single{si};
            try
                tic;
                tt = poly2qtt(funs_coef, d, 'method', method);
                single_runtime(di, rr, si) = toc;
                tts_single{si} = tt;
                single_ok(si) = true;
                single_rank(di, rr, si) = TTrankMean(tt);
                single_storage(di, rr, si) = mem(tt);
            catch ME
                % leave NaN
            end
        end

        % Batch error computation (ground truth evaluated once)
        ok_idx = find(single_ok);
        if ~isempty(ok_idx)
            errs = qtt_poly_error(tts_single(ok_idx), coef_single, 'cheb', [-1,1]);
            for k = 1:numel(ok_idx)
                single_error(di, rr, ok_idx(k)) = errs(k);
            end
        end

        %% ===== JOINT MULTI-POLYNOMIAL (M=M_multi) =====
        funs_coef_m = {coef_multi, 'cheb'};
        tts_joint = cell(n_joint, 1);
        joint_ok = false(n_joint, 1);
        for ji = 1:n_joint
            method = methods_joint{ji};
            try
                tic;
                tt = poly2qtt(funs_coef_m, d, 'method', method, 'order', 'last');
                joint_runtime(di, rr, ji) = toc;
                tts_joint{ji} = tt;
                joint_ok(ji) = true;
                joint_rank(di, rr, ji) = TTrankMean(tt);
                joint_storage(di, rr, ji) = mem(tt);
            catch ME
                % leave NaN
            end
        end

        % Batch error computation (ground truth evaluated once)
        ok_idx = find(joint_ok);
        if ~isempty(ok_idx)
            errs = qtt_poly_error(tts_joint(ok_idx), coef_multi, 'cheb', [-1,1], 'order', 'last');
            for k = 1:numel(ok_idx)
                joint_error(di, rr, ok_idx(k)) = errs(k);
            end
        end

        % Clear large temporaries to free memory
        clear tt tts_single tts_joint;

        % Progress indicator
        fprintf('  d=%2d done\n', d);
    end

    t_repeat = toc(t_repeat_start);
    t_total = toc(t_compute_start);
    fprintf('Repeat %d/%d finished in %.1f s (total elapsed: %.1f s / %.1f min)\n', ...
        rr, repeat, t_repeat, t_total, t_total/60);
end

%% Print per-d summary (mean +/- std across repeats)
for di = 1:n_d
    d = d_values(di);
    fprintf('\nSingle polynomial (M=1), d=%d, m=%d:\n', d, m);
    fprintf('%-20s %20s %18s %18s %20s\n', ...
        'Method', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
    fprintf('%s\n', repmat('-', 1, 98));
    for si = 1:n_single
        e = squeeze(single_error(di, :, si));
        t_r = squeeze(single_runtime(di, :, si));
        r = squeeze(single_rank(di, :, si));
        s = squeeze(single_storage(di, :, si));
        valid = ~isnan(e);
        if any(valid)
            fprintf('%-20s %9.2e ± %.2e %8d ± %5d %8.2f ± %5.2f %9.4f ± %.4f\n', ...
                labels_single{si}, mean(e(valid)), std(e(valid)), ...
                round(mean(s(valid))), round(std(s(valid))), ...
                mean(r(valid)), std(r(valid)), ...
                mean(t_r(valid)), std(t_r(valid)));
        else
            fprintf('%-20s %20s\n', labels_single{si}, 'ALL FAILED');
        end
    end

    fprintf('\nJoint polynomial (M=%d), d=%d, m=%d:\n', M_multi, d, m);
    fprintf('%-20s %20s %18s %18s %20s\n', ...
        'Method', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
    fprintf('%s\n', repmat('-', 1, 98));
    for ji = 1:n_joint
        e = squeeze(joint_error(di, :, ji));
        t_r = squeeze(joint_runtime(di, :, ji));
        r = squeeze(joint_rank(di, :, ji));
        s = squeeze(joint_storage(di, :, ji));
        valid = ~isnan(e);
        if any(valid)
            fprintf('%-20s %9.2e ± %.2e %8d ± %5d %8.2f ± %5.2f %9.4f ± %.4f\n', ...
                labels_joint{ji}, mean(e(valid)), std(e(valid)), ...
                round(mean(s(valid))), round(std(s(valid))), ...
                mean(r(valid)), std(r(valid)), ...
                mean(t_r(valid)), std(t_r(valid)));
        else
            fprintf('%-20s %20s\n', labels_joint{ji}, 'ALL FAILED');
        end
    end
end

% Save results
if ~exist('data', 'dir'), mkdir('data'); end
save(datafile, 'd_values', 'm', 'M_multi', 'repeat', ...
    'methods_single', 'labels_single', 'methods_joint', 'labels_joint', ...
    'single_error', 'single_runtime', 'single_rank', 'single_storage', ...
    'joint_error', 'joint_runtime', 'joint_rank', 'joint_storage', ...
    'n_single', 'n_joint', 'n_d');
fprintf('\nResults saved to: %s\n', datafile);

end % do_compute

%% Plotting & tables
if do_plot

if ~do_compute
    fprintf('Loading data from: %s\n', datafile);
    S = load(datafile);
    d_values = S.d_values;  m = S.m;  M_multi = S.M_multi;
    methods_single = S.methods_single;  labels_single = S.labels_single;
    methods_joint = S.methods_joint;    labels_joint = S.labels_joint;
    single_error = S.single_error;      single_runtime = S.single_runtime;
    single_rank = S.single_rank;        single_storage = S.single_storage;
    joint_error = S.joint_error;        joint_runtime = S.joint_runtime;
    joint_rank = S.joint_rank;          joint_storage = S.joint_storage;
    n_single = S.n_single;  n_joint = S.n_joint;
    n_d = S.n_d;  repeat = S.repeat;
end

figdir = 'figures/';
if saveresults && ~exist(figdir, 'dir')
    mkdir(figdir);
end

% Julia color scheme
julia_colors = [
    0.0, 0.6056, 0.9787;      % blue
    0.8889, 0.4356, 0.2781;   % orange
    0.2422, 0.6433, 0.3044;   % green
    0.7644, 0.4441, 0.8243;   % purple
];
% Single: LR=purple, FR=green, TTCross=blue, Horner=orange
single_colors = [julia_colors(4,:); julia_colors(3,:); julia_colors(1,:); julia_colors(2,:)];
% Joint: LR=purple, TTCross=blue
joint_colors  = [julia_colors(4,:); julia_colors(1,:)];

single_markers = {'o', 's', '^', 'd'};
joint_markers  = {'o', '^'};

%% --- Compute mean and std across repeats ---
% reshape ensures correct dimensions even when n_d==1
se_mean = reshape(mean(single_error, 2, 'omitnan'), [n_d, n_single]);
se_std  = reshape(std(single_error, 0, 2, 'omitnan'), [n_d, n_single]);
sr_mean = reshape(mean(single_rank, 2, 'omitnan'), [n_d, n_single]);
sr_std  = reshape(std(single_rank, 0, 2, 'omitnan'), [n_d, n_single]);
ss_mean = reshape(mean(single_storage, 2, 'omitnan'), [n_d, n_single]);
ss_std  = reshape(std(single_storage, 0, 2, 'omitnan'), [n_d, n_single]);
st_mean = reshape(mean(single_runtime, 2, 'omitnan'), [n_d, n_single]);
st_std  = reshape(std(single_runtime, 0, 2, 'omitnan'), [n_d, n_single]);

je_mean = reshape(mean(joint_error, 2, 'omitnan'), [n_d, n_joint]);
je_std  = reshape(std(joint_error, 0, 2, 'omitnan'), [n_d, n_joint]);
jr_mean = reshape(mean(joint_rank, 2, 'omitnan'), [n_d, n_joint]);
jr_std  = reshape(std(joint_rank, 0, 2, 'omitnan'), [n_d, n_joint]);
js_mean = reshape(mean(joint_storage, 2, 'omitnan'), [n_d, n_joint]);
js_std  = reshape(std(joint_storage, 0, 2, 'omitnan'), [n_d, n_joint]);
jt_mean = reshape(mean(joint_runtime, 2, 'omitnan'), [n_d, n_joint]);
jt_std  = reshape(std(joint_runtime, 0, 2, 'omitnan'), [n_d, n_joint]);

% Geometric mean and geometric std for error (symmetric on log scale)
% geo_mean = exp(mean(log(x))), geo_std_factor = exp(std(log(x)))
% error bars: [geo_mean/geo_std_factor, geo_mean*geo_std_factor]
se_log = log(single_error);
se_gmean = reshape(exp(mean(se_log, 2, 'omitnan')), [n_d, n_single]);
se_gstdf = reshape(exp(std(se_log, 0, 2, 'omitnan')), [n_d, n_single]);

je_log = log(joint_error);
je_gmean = reshape(exp(mean(je_log, 2, 'omitnan')), [n_d, n_joint]);
je_gstdf = reshape(exp(std(je_log, 0, 2, 'omitnan')), [n_d, n_joint]);

% Indices: all methods for error/runtime; skip FR (index 2) for rank/storage
plot_all_idx = 1:n_single;    % LR, FR, TTCross, Horner
plot_no_fr_idx = [1, 3, 4];   % LR, TTCross, Horner

%% --- Combined scaling figure (2 columns x 4 rows) ---
% Left column: single polynomial (N=1)
% Right column: joint polynomial (N=50)
% Rows: error, rank, storage, runtime
% Matches paper Figure 5 (scaling.tex)

fig = figure('Position', [100, 100, 800, 900]);

% ---- Row 1: Error (geometric mean with shaded bands) ----
% Single (left)
ax = subplot(4,2,1); hold on;
for si = plot_all_idx
    valid = ~isnan(se_gmean(:, si));
    if any(valid)
        gm = se_gmean(valid, si);
        gf = se_gstdf(valid, si);
        plot_shaded(ax, d_values(valid), gm, gm.*gf, gm./gf, ...
            single_colors(si,:), single_markers{si}, labels_single{si});
    end
end
hold off;
ylabel('$\ell_\infty$ error', 'Interpreter', 'latex');
title(sprintf('Single polynomial ($N = 1$)'), 'Interpreter', 'latex');
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([1e-13, 1e-6]);
set(gca, 'YTick', 10.^(-14:2:-6));
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;
legend('Location', 'best', 'FontSize', 7);

% Joint (right)
ax = subplot(4,2,2); hold on;
for ji = 1:n_joint
    valid = ~isnan(je_gmean(:, ji));
    if any(valid)
        gm = je_gmean(valid, ji);
        gf = je_gstdf(valid, ji);
        plot_shaded(ax, d_values(valid), gm, gm.*gf, gm./gf, ...
            joint_colors(ji,:), joint_markers{ji}, labels_joint{ji});
    end
end
hold off;
title(sprintf('Joint approximation ($N = %d$)', M_multi), 'Interpreter', 'latex');
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([1e-12, 1e3]);
set(gca, 'YTick', [1e-10, 1e-6, 1e-2, 1e2]);
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;
legend('Location', 'best', 'FontSize', 7);

% ---- Row 2: Rank (skip FR for single; all for joint) ----
% Single (left) — skip FR (index 2)
ax = subplot(4,2,3); hold on;
for si = plot_no_fr_idx
    valid = ~isnan(sr_mean(:, si));
    if any(valid)
        plot_shaded(ax, d_values(valid), sr_mean(valid, si), ...
            sr_mean(valid, si) + sr_std(valid, si), ...
            max(sr_mean(valid, si) - sr_std(valid, si), 0), ...
            single_colors(si,:), single_markers{si}, labels_single{si});
    end
end
hold off;
ylabel('avg. rank $\bar{r}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 9);
ylim([8, 18]);
set(gca, 'YTick', [10, 12, 14, 16]);
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;

% Joint (right)
ax = subplot(4,2,4); hold on;
for ji = 1:n_joint
    valid = ~isnan(jr_mean(:, ji));
    if any(valid)
        plot_shaded(ax, d_values(valid), jr_mean(valid, ji), ...
            jr_mean(valid, ji) + jr_std(valid, ji), ...
            max(jr_mean(valid, ji) - jr_std(valid, ji), 0), ...
            joint_colors(ji,:), joint_markers{ji}, labels_joint{ji});
    end
end
hold off;
set(gca, 'FontSize', 9);
ylim([5, 55]);
set(gca, 'YTick', [10, 20, 30, 40, 50]);
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;

% ---- Row 3: Storage (skip FR for single) ----
% Single (left)
ax = subplot(4,2,5); hold on;
for si = plot_no_fr_idx
    valid = ~isnan(ss_mean(:, si));
    if any(valid)
        plot_shaded(ax, d_values(valid), ss_mean(valid, si), ...
            ss_mean(valid, si) + ss_std(valid, si), ...
            max(ss_mean(valid, si) - ss_std(valid, si), 0), ...
            single_colors(si,:), single_markers{si}, labels_single{si});
    end
end
hold off;
ylabel('storage', 'Interpreter', 'latex');
set(gca, 'FontSize', 9);
ylim([-500, 10500]);
set(gca, 'YTick', [2000, 4000, 6000, 8000]);
set(gca, 'YTickLabel', {'2k', '4k', '6k', '8k'});
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;

% Joint (right)
ax = subplot(4,2,6); hold on;
for ji = 1:n_joint
    valid = ~isnan(js_mean(:, ji));
    if any(valid)
        plot_shaded(ax, d_values(valid), js_mean(valid, ji), ...
            js_mean(valid, ji) + js_std(valid, ji), ...
            max(js_mean(valid, ji) - js_std(valid, ji), 0), ...
            joint_colors(ji,:), joint_markers{ji}, labels_joint{ji});
    end
end
hold off;
set(gca, 'FontSize', 9);
ylim([0, 75000]);
set(gca, 'YTick', [15000, 30000, 45000, 60000]);
set(gca, 'YTickLabel', {'15k', '30k', '45k', '60k'});
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30, 'XTickLabel', []);
grid on; box on;

% ---- Row 4: Runtime (all methods) ----
% Single (left)
ax = subplot(4,2,7); hold on;
for si = plot_all_idx
    valid = ~isnan(st_mean(:, si));
    if any(valid)
        mu = st_mean(valid, si);
        sd = st_std(valid, si);
        plot_shaded(ax, d_values(valid), mu, mu + sd, max(mu - sd, eps), ...
            single_colors(si,:), single_markers{si}, labels_single{si});
    end
end
hold off;
xlabel('$d$ (grid size $2^d$)', 'Interpreter', 'latex');
ylabel('runtime (s)', 'Interpreter', 'latex');
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([1e-2, 2e0]);
set(gca, 'YTick', [1e-2, 1e-1, 1e0]);
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30);
grid on; box on;

% Joint (right)
ax = subplot(4,2,8); hold on;
for ji = 1:n_joint
    valid = ~isnan(jt_mean(:, ji));
    if any(valid)
        mu = jt_mean(valid, ji);
        sd = jt_std(valid, ji);
        plot_shaded(ax, d_values(valid), mu, mu + sd, max(mu - sd, eps), ...
            joint_colors(ji,:), joint_markers{ji}, labels_joint{ji});
    end
end
hold off;
xlabel('$d$ (grid size $2^d$)', 'Interpreter', 'latex');
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([4e-2, 2e0]);
set(gca, 'YTick', [1e-1, 1e0]);
xlim([d_values(1)-1, d_values(end)+1]);
set(gca, 'XTick', 10:4:30);
grid on; box on;

if saveresults
    filename = fullfile(figdir, 'scaling.pdf');
    save_figure_tight(fig, filename);
    fprintf('Saved: %s\n', filename);
end

%% --- Print summary tables ---
fprintf('\n=========================================================================\n');
fprintf('                    SCALING RESULTS SUMMARY                              \n');
fprintf('=========================================================================\n');
fprintf('  (mean ± std over %d repeats)\n', repeat);

fprintf('\n--- Single Polynomial (M=1), m=%d ---\n', m);
fprintf('%6s  %-20s %20s %18s %18s %20s\n', ...
    'd', 'Method', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
fprintf('%s\n', repmat('-', 1, 100));
for di = 1:n_d
    for si = 1:n_single
        if isnan(se_mean(di, si))
            fprintf('%6d  %-20s %20s\n', d_values(di), labels_single{si}, 'FAILED');
        else
            fprintf('%6d  %-20s %9.2e ± %.2e %8d ± %5d %8.2f ± %5.2f %9.4f ± %.4f\n', ...
                d_values(di), labels_single{si}, ...
                se_mean(di, si), se_std(di, si), ...
                round(ss_mean(di, si)), round(ss_std(di, si)), ...
                sr_mean(di, si), sr_std(di, si), ...
                st_mean(di, si), st_std(di, si));
        end
    end
end

fprintf('\n--- Joint Polynomial (M=%d), m=%d ---\n', M_multi, m);
fprintf('%6s  %-20s %20s %18s %18s %20s\n', ...
    'd', 'Method', 'L∞ Error', 'Storage', 'Avg Rank', 'Runtime (s)');
fprintf('%s\n', repmat('-', 1, 100));
for di = 1:n_d
    for ji = 1:n_joint
        if isnan(je_mean(di, ji))
            fprintf('%6d  %-20s %20s\n', d_values(di), labels_joint{ji}, 'FAILED');
        else
            fprintf('%6d  %-20s %9.2e ± %.2e %8d ± %5d %8.2f ± %5.2f %9.4f ± %.4f\n', ...
                d_values(di), labels_joint{ji}, ...
                je_mean(di, ji), je_std(di, ji), ...
                round(js_mean(di, ji)), round(js_std(di, ji)), ...
                jr_mean(di, ji), jr_std(di, ji), ...
                jt_mean(di, ji), jt_std(di, ji));
        end
    end
end



end % do_plot

t_wall = toc(t_wall_start);
fprintf('\n=== Scaling experiment complete (%.1f seconds / %.1f minutes) ===\n', t_wall, t_wall/60);
end

function plot_shaded(ax, x, y_center, y_upper, y_lower, color, marker, name)
%PLOT_SHADED Plot a line with a semi-transparent shaded band.
    x = x(:); y_center = y_center(:); y_upper = y_upper(:); y_lower = y_lower(:);
    x_fill = [x; flipud(x)];
    y_fill = [y_upper; flipud(y_lower)];
    fill(ax, x_fill, y_fill, color, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    plot(ax, x, y_center, ['-' marker], 'Color', color, 'LineWidth', 1.5, ...
         'MarkerSize', 5, 'MarkerFaceColor', color, 'DisplayName', name);
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
