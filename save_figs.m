function save_figs(folder)
% Saves every open figure as a PNG (light theme, white background) in "folder".
if ~exist(folder, 'dir'), mkdir(folder); end
fh = findobj('Type', 'figure');
for k = 1:numel(fh)
    try, theme(fh(k), 'light'); catch, end          % force light theme (R2025a+)
    name = regexprep(fh(k).Name, '[^\w]+', '_');
    if isempty(name), name = 'figure'; end
    exportgraphics(fh(k), fullfile(folder, sprintf('Fig%02d_%s.png', fh(k).Number, name)), ...
                   'Resolution', 200, 'BackgroundColor', 'white');
end
fprintf('Saved %d figures in %s\n', numel(fh), folder);
end
