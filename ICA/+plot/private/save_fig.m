function save_fig(fig, basePath, suffix)
if isempty(basePath), return; end
d = fileparts(basePath);
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end
fname = [basePath '_' suffix '.png'];
print(fig, fname, '-dpng', '-r150');
fprintf('ica.plot: saved %s\n', fname);
end
