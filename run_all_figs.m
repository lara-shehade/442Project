%% Runs every script, saves all figures (light mode) and zips them.
% The Command Window output of all scripts is also saved in run_log.txt
% (each script starts with clc, so only the last one stays on screen).
if exist('run_log.txt', 'file'), delete('run_log.txt'); end
diary('run_log.txt');
compare_1;      save_figs('figs_1_1');
fourierApp1_2;  save_figs('figs_1_2');
part2;          save_figs('figs_part2');      % plays sound for about 15 s
part3;          save_figs('figs_part3');
zip('all_figs.zip', {'figs_1_1', 'figs_1_2', 'figs_part2', 'figs_part3'});
disp('Done: all_figs.zip (figures) and run_log.txt (printed results)');
diary off
