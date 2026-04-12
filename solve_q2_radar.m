clear;
clc;
close all;

workdir = fileparts(mfilename('fullpath'));

radar_path = find_radar_mat(workdir);
if isempty(radar_path)
    summary = struct();
    summary.radar_dataset_found = false;
    summary.message = 'Place cazac_radar_dataset.mat (or similar name) in workspace root.';

    save(fullfile(workdir, 'q2_summary.mat'), 'summary');
    disp('Q2 summary:');
    disp(summary);
    return;
end

data = load(radar_path);
if ~isfield(data, 's') || ~isfield(data, 'r')
    error('MAT file must contain variables s and r.');
end

s = squeeze(data.s);
r = squeeze(data.r);
s = complex(double(s(:)));
r = complex(double(r(:)));

n_clutter = 500;
v_train = r(1:n_clutter);
alpha_hat = estimate_ar1_alpha(v_train);

% FIR pre-whitener: H_w(z) = 1 - alpha z^-1
b = [1, -alpha_hat];
a = 1;

% Receiver A: naive matched filter / cross-correlation
[c_a, lags_a, c_a_db] = xcorr_db(r, s);
[~, peak_idx_a] = max(abs(c_a));
peak_lag_a = lags_a(peak_idx_a);
snr_a_db = output_snr_db(c_a, peak_idx_a, 50);

% Receiver B: pre-whiten both then correlate
r_w = filter(b, a, r);
s_w = filter(b, a, s);
[c_b, lags_b, c_b_db] = xcorr_db(r_w, s_w);
[~, peak_idx_b] = max(abs(c_b));
peak_lag_b = lags_b(peak_idx_b);
snr_b_db = output_snr_db(c_b, peak_idx_b, 50);

% PSD evidence for whitening performance on clutter-only segment
v_w = filter(b, a, v_train);
[f_psd, psd_v_db] = simple_psd_db(v_train, 2048);
[~, psd_vw_db] = simple_psd_db(v_w, 2048);

f1 = figure('Visible', 'off');
plot(lags_a, c_a_db, 'LineWidth', 1.1); hold on;
xline(peak_lag_a, '--k', sprintf('Peak lag=%d', peak_lag_a));
grid on;
ylim([-100, 2]);
xlabel('Lag (samples)');
ylabel('|Cross-correlation| (dB, normalized)');
title('Q2 Receiver A: Raw cross-correlation');
legend({'Receiver A', 'Peak lag'}, 'Location', 'best');
saveas(f1, fullfile(workdir, 'q2_receiver_A_xcorr_db.png'));
close(f1);

f2 = figure('Visible', 'off');
plot(lags_b, c_b_db, 'LineWidth', 1.1); hold on;
xline(peak_lag_b, '--k', sprintf('Peak lag=%d', peak_lag_b));
grid on;
ylim([-100, 2]);
xlabel('Lag (samples)');
ylabel('|Cross-correlation| (dB, normalized)');
title('Q2 Receiver B: Whitened matched filter');
legend({'Receiver B', 'Peak lag'}, 'Location', 'best');
saveas(f2, fullfile(workdir, 'q2_receiver_B_xcorr_db.png'));
close(f2);

f3 = figure('Visible', 'off');
plot(f_psd, psd_v_db, 'LineWidth', 1.2); hold on;
plot(f_psd, psd_vw_db, 'LineWidth', 1.2);
grid on;
xlabel('Normalized frequency (cycles/sample)');
ylabel('PSD (dB, normalized)');
title('Q2: Clutter PSD before and after pre-whitening');
legend({'Raw clutter v[n]', 'Whitened clutter v_w[n]'}, 'Location', 'best');
saveas(f3, fullfile(workdir, 'q2_psd_before_after_whitening.png'));
close(f3);

f4 = figure('Visible', 'off');
plot(lags_a, c_a_db, 'LineWidth', 1.0); hold on;
plot(lags_b, c_b_db, 'LineWidth', 1.0);
xline(peak_lag_a, '--k', sprintf('Peak lag=%d', peak_lag_a));
grid on;
ylim([-100, 2]);
xlabel('Lag (samples)');
ylabel('|Cross-correlation| (dB, normalized)');
title('Q2: Receiver A vs Receiver B comparison');
legend({'Receiver A (raw)', 'Receiver B (whitened)', 'Detected lag'}, 'Location', 'best');
saveas(f4, fullfile(workdir, 'q2_receivers_overlay_db.png'));
close(f4);

summary = struct();
summary.radar_dataset_found = true;
summary.radar_dataset_path = radar_path;
summary.alpha_hat_real = real(alpha_hat);
summary.alpha_hat_imag = imag(alpha_hat);
summary.receiver_A_peak_lag = peak_lag_a;
summary.receiver_A_output_snr_db = snr_a_db;
summary.receiver_B_peak_lag = peak_lag_b;
summary.receiver_B_output_snr_db = snr_b_db;
summary.snr_gain_db_B_minus_A = snr_b_db - snr_a_db;

save(fullfile(workdir, 'q2_summary.mat'), 'summary');
disp('Q2 summary:');
disp(summary);


function radar_file = find_radar_mat(workdir)
    radar_file = '';
    mats = dir(fullfile(workdir, '*.mat'));
    for i = 1:numel(mats)
        name = lower(strrep(mats(i).name, '_', ' '));
        if contains(name, 'cazac') && contains(name, 'radar') && contains(name, 'dataset')
            radar_file = mats(i).name;
            return;
        end
    end
end

function alpha_hat = estimate_ar1_alpha(v)
    num = sum(v(2:end) .* conj(v(1:end-1)));
    den = sum(abs(v(1:end-1)).^2);
    alpha_hat = num / den;
end

function [c, lags, mag_db] = xcorr_db(x, y)
    [c, lags] = xcorr(x, y);
    mag = abs(c);
    mag_db = 20 * log10(mag ./ (max(mag) + eps) + eps);
end

function snr_db = output_snr_db(corr_vec, peak_idx, exclusion)
    mag2 = abs(corr_vec).^2;
    signal_power = mag2(peak_idx);

    mask = true(size(mag2));
    lo = max(1, peak_idx - exclusion);
    hi = min(numel(mask), peak_idx + exclusion);
    mask(lo:hi) = false;

    noise_power = mean(mag2(mask));
    snr_db = 10 * log10(signal_power / noise_power);
end

function [f, p_db] = simple_psd_db(x, nfft)
    x = x(:);
    x = x - mean(x);
    X = fftshift(fft(x, nfft));
    p = abs(X).^2;
    p_db = 10 * log10(p ./ (max(p) + eps) + eps);
    f = linspace(-0.5, 0.5, nfft);
end
