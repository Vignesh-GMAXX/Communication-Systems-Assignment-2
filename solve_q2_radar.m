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

% Absolute dB views used by the additional presentation plots
c_a_abs_db = 20 * log10(abs(c_a) + eps);
c_b_abs_db = 20 * log10(abs(c_b) + eps);

% PSD evidence for whitening performance on clutter-only segment
v_w = filter(b, a, v_train);
[f_psd, psd_v_db] = simple_psd_db(v_train, 2048);
[~, psd_vw_db] = simple_psd_db(v_w, 2048);

nfft_rad = 1024;
omega = linspace(-pi, pi, nfft_rad);
v_fft = fftshift(fft(v_train, nfft_rad));
vw_fft = fftshift(fft(v_w, nfft_rad));

psd_raw_db = 10 * log10(abs(v_fft).^2 / nfft_rad + eps);
psd_white_db = 10 * log10(abs(vw_fft).^2 / nfft_rad + eps);

% Set reference to whitened mean level so whitened PSD appears around 0 dB.
psd_ref_db = mean(psd_white_db);
psd_raw_rel_db = psd_raw_db - psd_ref_db;
psd_white_rel_db = psd_white_db - psd_ref_db;

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

% Additional plot 1: show only lag window [1, 1000] in absolute dB
lag_window_mask_a = lags_a >= 1 & lags_a <= 1000;
lag_window_mask_b = lags_b >= 1 & lags_b <= 1000;

lag_window_a = lags_a(lag_window_mask_a);
lag_window_b = lags_b(lag_window_mask_b);
c_a_window_db = c_a_abs_db(lag_window_mask_a);
c_b_window_db = c_b_abs_db(lag_window_mask_b);

window_ymin = max(0, floor(min([c_a_window_db; c_b_window_db])) - 2);
window_ymax = ceil(max([c_a_window_db; c_b_window_db])) + 2;

f5 = figure('Visible', 'off');
subplot(2, 1, 1);
plot(lag_window_a, c_a_window_db, 'b', 'LineWidth', 1.0);
grid on;
xlim([1, 1000]);
ylim([window_ymin, window_ymax]);
xlabel('Lag');
ylabel('Magnitude (dB)');
title('Receiver A: Naive Matched Filter Output');

subplot(2, 1, 2);
plot(lag_window_b, c_b_window_db, 'r', 'LineWidth', 1.0);
grid on;
xlim([1, 1000]);
ylim([window_ymin, window_ymax]);
xlabel('Lag');
ylabel('Magnitude (dB)');
title('Receiver B: Optimal Whitened Matched Filter Output');

saveas(f5, fullfile(workdir, 'q2_receivers_window_1_1000_db.png'));
close(f5);

% Additional plot 2: PSD in normalized rad/sample (x in [-pi, pi])
f6 = figure('Visible', 'off');
plot(omega, psd_raw_rel_db, 'b', 'LineWidth', 1.0); hold on;
plot(omega, psd_white_rel_db, 'r', 'LineWidth', 1.0);
grid on;
xlim([-pi, pi]);
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Power/Frequency (dB/(rad/sample))');
title('Power Spectral Density: Raw vs. Whitened Signal');
legend({'Raw Signal r[n] (Colored Clutter)', 'Whitened Signal r_w[n] (Flattened)'}, 'Location', 'northeast');
saveas(f6, fullfile(workdir, 'q2_psd_raw_vs_whitened_rad.png'));
close(f6);

% Additional plot 3: extended lag axis with zero padding outside valid lags
lag_ext = -5000:5000;
c_a_ext = zeros(size(lag_ext));
c_b_ext = zeros(size(lag_ext));

[is_lag_a, idx_lag_a] = ismember(lags_a, lag_ext);
[is_lag_b, idx_lag_b] = ismember(lags_b, lag_ext);
c_a_ext(idx_lag_a(is_lag_a)) = c_a(is_lag_a);
c_b_ext(idx_lag_b(is_lag_b)) = c_b(is_lag_b);

floor_mag = 10^(-140/20);
c_a_ext_db = 20 * log10(abs(c_a_ext) + floor_mag);
c_b_ext_db = 20 * log10(abs(c_b_ext) + floor_mag);

f7 = figure('Visible', 'off');
subplot(2, 1, 1);
plot(lag_ext, c_a_ext_db, 'b', 'LineWidth', 1.0);
grid on;
xlim([-5000, 5000]);
ylim([-150, 50]);
xlabel('Lag');
ylabel('Magnitude (dB)');
title('Receiver A: Naive Matched Filter Output');

subplot(2, 1, 2);
plot(lag_ext, c_b_ext_db, 'r', 'LineWidth', 1.0);
grid on;
xlim([-5000, 5000]);
ylim([-150, 50]);
xlabel('Lag');
ylabel('Magnitude (dB)');
title('Receiver B: Optimal Whitened Matched Filter Output');

saveas(f7, fullfile(workdir, 'q2_receivers_extended_lag_db.png'));
close(f7);

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
summary.generated_plot_files = {
    'q2_receiver_A_xcorr_db.png', ...
    'q2_receiver_B_xcorr_db.png', ...
    'q2_psd_before_after_whitening.png', ...
    'q2_receivers_overlay_db.png', ...
    'q2_receivers_window_1_1000_db.png', ...
    'q2_psd_raw_vs_whitened_rad.png', ...
    'q2_receivers_extended_lag_db.png'
};

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
