function save_fig(fig, basePath, suffix)
if isempty(basePath), return; end
d = fileparts(basePath);
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end
[~, annotLabel] = fileparts(basePath);
annotation('textbox', [0.01 0.01, 0.01 0.01], 'String', basePath, ...
    'FontSize', 6, 'LineStyle', 'none', 'Interpreter', 'none');
fname = [basePath '_' suffix '.png'];
if gedai.printFigure(fig, fname)
    fprintf('gedai.eval: saved %s\n', fname);
end
end
