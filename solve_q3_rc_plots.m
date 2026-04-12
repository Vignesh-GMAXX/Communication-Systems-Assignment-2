clear;
clc;
close all;

workdir = fileparts(mfilename('fullpath'));

T = 1;
beta_list = [0, 0.25, 0.5, 1.0];

t = linspace(-8*T, 8*T, 4001);
f = linspace(-1.5, 1.5, 4001);

% Time-domain RC impulse response comparison
f1 = figure('Visible', 'off');
for k = 1:numel(beta_list)
    b = beta_list(k);
    h = rc_impulse(t, T, b);
    plot(t/T, h, 'LineWidth', 1.2); hold on;
end
grid on;
xlabel('t/T');
ylabel('h_{RC}(t)');
title('Q3: Raised-cosine impulse response for different \beta');
legend({'\beta=0', '\beta=0.25', '\beta=0.5', '\beta=1.0'}, 'Location', 'best');
saveas(f1, fullfile(workdir, 'q3_rc_impulse_beta_comparison.png'));
close(f1);

% Frequency-domain RC response comparison
f2 = figure('Visible', 'off');
for k = 1:numel(beta_list)
    b = beta_list(k);
    H = rc_frequency(f, T, b);
    plot(f, H, 'LineWidth', 1.2); hold on;
end
grid on;
xlabel('f (normalized, T=1)');
ylabel('H_{RC}(f)');
title('Q3: Raised-cosine spectrum for different \beta');
legend({'\beta=0', '\beta=0.25', '\beta=0.5', '\beta=1.0'}, 'Location', 'best');
axis([-1.2 1.2 -0.05 1.05]);
saveas(f2, fullfile(workdir, 'q3_rc_frequency_beta_comparison.png'));
close(f2);

fprintf('Generated Q3 figures:\n');
fprintf(' - q3_rc_impulse_beta_comparison.png\n');
fprintf(' - q3_rc_frequency_beta_comparison.png\n');


function h = rc_impulse(t, T, beta)
    x = t ./ T;

    if beta == 0
        h = sinc_local(x);
        return;
    end

    num = sinc_local(x) .* cos(pi * beta * x);
    den = 1 - (2 * beta * x).^2;
    h = num ./ den;

    % Handle removable singularities at x = +-1/(2*beta)
    idx = abs(abs(x) - 1/(2*beta)) < 1e-10;
    if any(idx)
        h(idx) = (beta/2) * sin(pi/(2*beta));
    end

    % Handle x = 0 explicitly for numerical stability
    idx0 = abs(x) < 1e-12;
    h(idx0) = 1;
end

function y = sinc_local(x)
    y = ones(size(x));
    nz = abs(x) > 1e-12;
    y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end

function H = rc_frequency(f, T, beta)
    af = abs(f);
    H = zeros(size(f));

    if beta == 0
        H(af <= 1/(2*T)) = T;
        return;
    end

    f1 = (1-beta)/(2*T);
    f2 = (1+beta)/(2*T);

    H(af <= f1) = T;

    idx = (af > f1) & (af <= f2);
    H(idx) = (T/2) .* (1 + cos((pi*T/beta) .* (af(idx) - f1)));
end
