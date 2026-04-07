clear;
clc;
close all;

workdir = fileparts(mfilename('fullpath'));

bits = 8;
gamma = linspace(0.5, 7.0, 500);

[p_gran, p_ol, p_tot] = theoretical_noise_curves(gamma, bits, 1.0);
[~, idx_opt] = min(p_tot);
gamma_opt = gamma(idx_opt);
sqnr_theory_db = 10*log10(1./p_tot);

f1 = figure('Visible', 'off');
plot(gamma, p_gran, 'LineWidth', 1.0); hold on;
plot(gamma, p_ol, 'LineWidth', 1.0);
plot(gamma, p_tot, 'LineWidth', 2.0);
xline(gamma_opt, '--k', sprintf('gamma_{opt}=%.3f', gamma_opt));
set(gca, 'YScale', 'log');
grid on;
xlabel('Loading factor gamma = Vmax/sigma');
ylabel('Noise power');
title('Q1: Theoretical quantization-noise components');
legend({'Theoretical granular noise', 'Theoretical overload noise', 'Theoretical total noise', 'gamma opt'}, 'Location', 'best');
saveas(f1, fullfile(workdir, 'q1_theoretical_noise_vs_gamma.png'));
close(f1);

f2 = figure('Visible', 'off');
plot(gamma, sqnr_theory_db, 'LineWidth', 2.0); hold on;
xline(gamma_opt, '--k', sprintf('gamma_{opt}=%.3f', gamma_opt));
grid on;
xlabel('Loading factor gamma = Vmax/sigma');
ylabel('Theoretical SQNR (dB)');
title('Q1: Theoretical SQNR vs loading factor');
legend({'Theoretical SQNR', 'gamma opt'}, 'Location', 'best');
saveas(f2, fullfile(workdir, 'q1_theoretical_sqnr_vs_gamma.png'));
close(f2);

q1_summary = struct();
q1_summary.bits = bits;
q1_summary.gamma_opt_theory = gamma_opt;
q1_summary.sqnr_opt_theory_db = sqnr_theory_db(idx_opt);
q1_summary.ofdm_dataset_found = false;

ofdm_path = find_ofdm_mat(workdir);
if ~isempty(ofdm_path)
    x = load_first_vector_from_mat(ofdm_path);
    x = x - mean(x);
    sigma_emp = std(x);

    sqnr_emp_db = zeros(size(gamma));
    for i = 1:numel(gamma)
        vmax = gamma(i) * sigma_emp;
        xq = quantize_midtread(x, vmax, bits);
        sqnr_emp_db(i) = empirical_sqnr(x, xq);
    end

    [~, idx_emp] = max(sqnr_emp_db);

    f3 = figure('Visible', 'off');
    plot(gamma, sqnr_theory_db, 'LineWidth', 1.5); hold on;
    plot(gamma, sqnr_emp_db, 'LineWidth', 1.5);
    xline(gamma_opt, '--k', sprintf('Theory gamma_{opt}=%.3f', gamma_opt));
    xline(gamma(idx_emp), ':r', sprintf('Empirical gamma_{opt}=%.3f', gamma(idx_emp)));
    grid on;
    xlabel('Loading factor gamma = Vmax/sigma');
    ylabel('SQNR (dB)');
    title('Q1: Theoretical vs empirical SQNR');
    legend({'Theoretical SQNR', 'Empirical SQNR (OFDM dataset)', 'Theory gamma opt', 'Empirical gamma opt'}, 'Location', 'best');
    saveas(f3, fullfile(workdir, 'q1_theory_vs_empirical_sqnr.png'));
    close(f3);

    q1_summary.ofdm_dataset_found = true;
    q1_summary.ofdm_dataset_path = ofdm_path;
    q1_summary.sigma_empirical = sigma_emp;
    q1_summary.gamma_opt_empirical = gamma(idx_emp);
    q1_summary.sqnr_opt_empirical_db = sqnr_emp_db(idx_emp);
end

save(fullfile(workdir, 'q1_summary.mat'), 'q1_summary');
disp('Q1 summary:');
disp(q1_summary);


function y = qfunc_local(x)
    y = 0.5 * erfc(x ./ sqrt(2));
end

function y = phi(x)
    y = exp(-0.5 .* x.^2) ./ sqrt(2*pi);
end

function xq = quantize_midtread(x, vmax, bits)
    levels = 2^bits;
    delta = 2*vmax/(levels - 1);
    x_clip = min(max(x, -vmax), vmax);
    xq = delta .* round(x_clip ./ delta);
end

function [p_gran, p_ol, p_tot] = theoretical_noise_curves(gamma, bits, sigma)
    levels = 2^bits;
    vmax = gamma .* sigma;
    delta = 2 .* vmax ./ (levels - 1);

    p_in = 1 - 2 .* qfunc_local(gamma);
    p_gran = (delta.^2 / 12) .* p_in;

    p_ol = 2 * sigma^2 .* ((1 + gamma.^2) .* qfunc_local(gamma) - gamma .* phi(gamma));
    p_tot = p_gran + p_ol;
end

function ofdm_file = find_ofdm_mat(workdir)
    ofdm_file = '';
    mats = dir(fullfile(workdir, '*.mat'));
    for i = 1:numel(mats)
        name = lower(strrep(mats(i).name, '_', ' '));
        if contains(name, 'ofdm') && contains(name, 'dataset')
            ofdm_file = mats(i).name;
            return;
        end
    end
end

function x = load_first_vector_from_mat(filename)
    data = load(filename);
    f = fieldnames(data);
    x = [];
    for i = 1:numel(f)
        v = data.(f{i});
        if isnumeric(v)
            if isvector(v) && ~isempty(v)
                x = real(double(v(:)));
                return;
            end
            if ismatrix(v) && any(size(v) == 1)
                x = real(double(v(:)));
                return;
            end
        end
    end
    error('No usable 1-D numeric array found in %s', filename);
end

function s = empirical_sqnr(x, xq)
    signal_power = mean(abs(x).^2);
    noise_power = mean(abs(x - xq).^2);
    s = 10*log10(signal_power / noise_power);
end
