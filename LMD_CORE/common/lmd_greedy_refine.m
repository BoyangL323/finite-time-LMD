function [edges, info] = lmd_greedy_refine(wgrid, Jgrid, Nb, tgrid, wfun, n0, seed, rank_mode, i_lo, i_hi, err_mode, cfg)
%LMD_GREEDY_REFINE  E_b-driven adaptive bin refinement.
%
%   [edges, info] = lmd_greedy_refine(wgrid, Jgrid, Nb, tgrid, wfun, n0, ...
%                                     seed, rank_mode, i_lo, i_hi, err_mode)
%
% err_mode is passed straight to lmd_bin_error_exact and MUST match the seed:
% 'l2' for route A / legacy, 'phase' for route B, and 'phase_sys' for LMD-S.
% Refining a partition against another route's functional optimizes the wrong
% target and is prohibited by the caller.
%
% CAVEAT that applies to both functionals: E_b is a per-bin surrogate and the
% greedy loop minimizes sum_b E_b, whereas the true error is the norm of the
% COHERENT sum of the per-bin defects.  Cross-bin cancellation is therefore not
% modelled.  This is the most likely reason the greedy rule still loses to a
% plain equal-mass partition in some large-budget scans; run_article_ablation
% records this comparison in the current environment.  This behavior is
% it is a theory gap, not a tuning problem - do not paper over it with a knob.
%
%   1. start from n0 equal-mass bins on the chosen seed measure
%   2. compute Eb for every bin (cached by lmd_bin_error_exact)
%   3. split the selected bin at its M0-median
%   4. repeat until numel(edges)-1 == Nb
%
% rank_mode selects step 3:
%
%   'abs'     split the bin with the largest absolute Eb.  This is the
%             historically unstable ranking: it piles modes into structural
%             peaks and can starve the high-frequency background.  It is kept
%             only for the versioned ablation, which writes its own metrics.
%
%   'permode' DEFAULT.  Split the bin with the largest ACTUAL error reduction
%             per added mode,
%                 gain_b = Eb - ( Eb_left + Eb_right ) ,
%             which costs exactly one extra mode.  A peak bin that is already
%             resolved has a small gain even though its absolute Eb is large, so
%             the background is no longer starved.  This is the fix required by
%             spec 3.5 problem 1 and it introduces NO new tunable parameter
%             (constraint C-2).
%
% There is deliberately no peak detection and no core/tail architecture
% (constraints C-4, C-5): structural peaks must be resolved by the error
% indicator alone.
%
% n0 defaults to min(Nb, 4).  It is an implementation constant, not a knob: the
% greedy loop refines down to the same budget regardless.
%
% info fields: Eb_final (per bin), Eb_total, gains, n_split, rank_mode,
%              i_edges (index-space boundaries), starved (true if the rule ran
%              out of splittable bins before reaching Nb).

wgrid = wgrid(:);
Jgrid = Jgrid(:);
if nargin < 5 || isempty(wfun),      wfun = @(t) ones(size(t)); end
if nargin < 6 || isempty(n0),        n0 = min(Nb, 4);           end
if nargin < 7 || isempty(seed),      seed = Jgrid;              end
if nargin < 8 || isempty(rank_mode), rank_mode = 'permode';     end
if nargin < 9 || isempty(i_lo),      i_lo = 1;                  end
if nargin < 10 || isempty(i_hi),     i_hi = numel(wgrid);       end
if nargin < 11 || isempty(err_mode), err_mode = 'l2';           end
if nargin < 12, cfg = []; end

rank_mode = lower(char(rank_mode));
Nb = max(1, round(Nb));
n0 = max(1, min(round(n0), Nb));

% ---- 1. initial equal-mass bins on the seed measure -------------------------
iq = local_equal_mass_indices(wgrid, seed, i_lo, i_hi, n0);

% ---- 2..4. greedy loop ------------------------------------------------------
starved = false;
while numel(iq) - 1 < Nb
    nb = numel(iq) - 1;
    score = -Inf(nb, 1);
    isplit = zeros(nb, 1);
    for b = 1:nb
        [score(b), isplit(b)] = local_bin_score(wgrid, Jgrid, iq(b), iq(b+1), ...
                                                tgrid, wfun, rank_mode, err_mode, cfg);
    end
    [best, ib] = max(score);
    if ~isfinite(best) || isplit(ib) <= iq(ib) || isplit(ib) >= iq(ib+1)
        starved = true;
        break;
    end
    iq = [iq(1:ib); isplit(ib); iq(ib+1:end)];
end

nb = numel(iq) - 1;
Eb_final = zeros(nb, 1);
for b = 1:nb
    Eb_final(b) = lmd_bin_error_exact(wgrid, Jgrid, iq(b), iq(b+1), tgrid, wfun, err_mode, cfg);
end

edges = wgrid(iq);
edges = edges(:);
edges(1)   = wgrid(i_lo);
edges(end) = wgrid(i_hi);

info = struct();
info.Eb_final  = Eb_final;
info.Eb_total  = sum(Eb_final);
info.n_bins    = nb;
info.n0        = n0;
info.rank_mode = rank_mode;
info.i_edges   = iq(:);
info.starved   = starved;
end

% ---------------------------------------------------------------------------
function [score, isplit] = local_bin_score(wgrid, Jgrid, i0, i1, tgrid, wfun, rank_mode, err_mode, cfg)
score  = -Inf;
isplit = 0;
if i1 - i0 < 2
    return;                                     % not splittable
end
Eb = lmd_bin_error_exact(wgrid, Jgrid, i0, i1, tgrid, wfun, err_mode, cfg);
im = local_mass_median_index(wgrid, Jgrid, i0, i1);
if im <= i0 || im >= i1
    return;
end
isplit = im;
switch rank_mode
    case 'abs'
        score = Eb;
    case 'permode'
        El = lmd_bin_error_exact(wgrid, Jgrid, i0, im, tgrid, wfun, err_mode, cfg);
        Er = lmd_bin_error_exact(wgrid, Jgrid, im, i1, tgrid, wfun, err_mode, cfg);
        score = Eb - (El + Er);                 % exactly one extra mode
    otherwise
        error('lmd_greedy_refine:BadRank', ...
              'cfg.lmd_greedy_rank must be ''abs'' or ''permode'', got %s.', rank_mode);
end
end

% ---------------------------------------------------------------------------
function im = local_mass_median_index(wgrid, Jgrid, i0, i1)
ww = wgrid(i0:i1);
JJ = Jgrid(i0:i1);
cum = cumtrapz(ww, JJ);
tot = cum(end);
if ~(tot > 0)
    im = i0 + floor((i1 - i0) / 2);
    return;
end
k = find(cum >= 0.5 * tot, 1, 'first');
if isempty(k), k = numel(ww); end
im = i0 + k - 1;
im = min(max(im, i0 + 1), i1 - 1);
end

% ---------------------------------------------------------------------------
function iq = local_equal_mass_indices(wgrid, dens, i_lo, i_hi, B)
ww = wgrid(i_lo:i_hi);
dd = max(dens(i_lo:i_hi), 0);
cum = cumtrapz(ww, dd);
tot = cum(end);
if ~(tot > 0) || B < 2
    iq = unique(round(linspace(i_lo, i_hi, max(B, 1) + 1))).';
    iq = iq(:);
    return;
end
targets = linspace(0, tot, B + 1).';
iq = zeros(B + 1, 1);
iq(1)   = i_lo;
iq(end) = i_hi;
for b = 2:B
    k = find(cum >= targets(b), 1, 'first');
    if isempty(k), k = numel(ww); end
    iq(b) = i_lo + k - 1;
end
% enforce strict monotonicity in index space
for b = 2:numel(iq)
    if iq(b) <= iq(b-1)
        iq(b) = iq(b-1) + 1;
    end
end
iq = min(iq, i_hi);
iq = unique(iq(:), 'stable');
if iq(end) ~= i_hi, iq(end+1) = i_hi; end
end
