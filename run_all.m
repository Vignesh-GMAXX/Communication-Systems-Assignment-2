clear;
clc;
close all;

plot_files = {
	'q1_theoretical_noise_vs_gamma.png', ...
	'q1_theoretical_sqnr_vs_gamma.png', ...
	'q1_theory_vs_empirical_sqnr.png', ...
	'q2_receiver_A_xcorr_db.png', ...
	'q2_receiver_B_xcorr_db.png'
};

for k = 1:numel(plot_files)
	if exist(plot_files{k}, 'file') == 2
		delete(plot_files{k});
	end
end

fprintf('\n=== Running solve_q1_quantizer.m ===\n');
run('solve_q1_quantizer.m');

fprintf('\n=== Running solve_q2_radar.m ===\n');
run('solve_q2_radar.m');

fprintf('\nAll MATLAB scripts finished.\n');
