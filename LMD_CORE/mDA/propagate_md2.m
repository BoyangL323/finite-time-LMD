function trajectory = propagate_md2(omega, lambda, cfg)
%PROPAGATE_MD2 Thin driver for the original DavydovD2Rabi TDVP function.
% The input lambda is the coupling amplitude. Time and frequency use reciprocal units.
arguments
    omega (:,1) double {mustBeReal, mustBeFinite}
    lambda (:,1) double {mustBeReal, mustBeFinite}
    cfg (1,1) struct
end
if isempty(omega) || numel(omega) ~= numel(lambda)
    error('propagate_md2:Modes', 'Provide one coupling per oscillator.');
end
modeCount = numel(omega);
multiplicity = cfg.MM;
systemSize = cfg.Ns;
hamiltonian = (cfg.epsilon * cfg.sigmaz - cfg.Delta * cfg.sigmax) / 2;

originalRng = rng;
restoreRng = onCleanup(@() rng(originalRng));
rng(cfg.rng_seed, 'twister');
spin = cfg.initial_system_state(:);
spin = spin / norm(spin);
amplitudes = zeros(multiplicity, systemSize);
amplitudes(1,:) = spin.';
amplitudes = amplitudes + cfg.init_A_noise * (2 * rand(size(amplitudes)) - 1);
displacements = cfg.init_F_noise * (2 * rand(multiplicity, modeCount) - 1);
state = [amplitudes(:); reshape(displacements.', [], 1)];
[~, initialNorm] = observabled2rwa(state, multiplicity, systemSize, modeCount, cfg.sigmaz);
amplitudes = amplitudes / sqrt(initialNorm);
state(1:numel(amplitudes)) = amplitudes(:);

times = linspace(0, cfg.t_end, cfg.n_t);
options = odeset('RelTol', cfg.RelTol, 'AbsTol', cfg.AbsTol);
if ~isempty(cfg.MaxStep)
    options = odeset(options, 'MaxStep', cfg.MaxStep);
end
[time, states] = ode113(@(t, y) DavydovD2Rabi(t, y, multiplicity, systemSize, ...
    modeCount, hamiltonian, cfg.VV1, lambda.', omega.'), times, state, options);

Pz = zeros(size(time));
normSquared = zeros(size(time));
for k = 1:numel(time)
    [Pz(k), normSquared(k)] = observabled2rwa( ...
        states(k,:).', multiplicity, systemSize, modeCount, cfg.sigmaz);
end
trajectory = table(time, (Pz + 1) / 2, normSquared, ...
    'VariableNames', {'time', 'P1', 'state_norm_squared'});
end
