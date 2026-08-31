function [sel, ok] = pickerDialog(itemLabels, preselected, dlgName, promptStr, confirmStr, classNames, classOfItem)
% PICKERDIALOG  Modal multi-select picker used by viewer.eegCompareViewer for
% all three of its list interactions: which channels to plot, which
% components to plot, and which components to subtract.
%
%   [sel, ok] = pickerDialog(itemLabels, preselected, dlgName, promptStr, ...
%                            confirmStr, classNames, classOfItem)
%
%   itemLabels  cellstr, one entry per selectable item
%   preselected indices highlighted when the dialog opens
%   dlgName     window title
%   promptStr   text above the list
%   confirmStr  label on the confirm button (e.g. 'Subtract')
%   classNames  cellstr of group names (IC classes, channel montages, ...),
%               {} to hide the group row
%   classOfItem nItems x nGroups logical membership matrix, or a numeric
%               per-item group index (0 = uncategorised) for exclusive groups
%
%   sel  selected indices (equal to preselected if cancelled)
%   ok   true if the user confirmed
%
% This is hand-rolled rather than a listdlg call because listdlg offers no
% category-based bulk selection ("highlight every Heart component") and no
% custom confirm label. Built from the same plain figure/uicontrol stack as
% the viewer itself -- deliberately not uifigure, see the renderer note in
% eegCompareViewer.m.
if nargin < 6, classNames  = {}; end
if nargin < 7, classOfItem = []; end

sel = preselected(:).';
ok  = false;

nItems     = numel(itemLabels);
sel        = sel(sel >= 1 & sel <= nItems);
hasClasses = ~isempty(classNames) && ~isempty(classOfItem);

% Group membership is a logical nItems x nGroups MATRIX, because groups are
% not mutually exclusive: an IC is both "Eye" and "bad" at the same time and
% either must be selectable. A numeric one-group-per-item vector is still
% accepted and expanded, for callers with a genuinely exclusive grouping.
if hasClasses
    if islogical(classOfItem)
        groupMask = reshape(classOfItem, nItems, []);
    else
        groupMask = false(nItems, numel(classNames));
        for gi = 1:numel(classNames)
            groupMask(:, gi) = classOfItem(:) == gi;
        end
    end
    hasClasses = size(groupMask, 2) == numel(classNames) && any(groupMask(:));
end

W = 380; H = 560;
margin = 14; inner = W - 2*margin;
if hasClasses, listBottom = 116; else, listBottom = 84; end

dlg = figure('Name', dlgName, 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'Color', 'w', 'Units', 'pixels', 'Resize', 'off', ...
    'Position', [100 100 W H], 'WindowStyle', 'modal', ...
    'CloseRequestFcn', @(s,e) onCancel());
movegui(dlg, 'center');

uicontrol(dlg, 'Style', 'text', 'String', promptStr, 'Units', 'pixels', ...
    'Position', [margin H-52 inner 38], 'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'w');

% Min/Max spread of more than 1 is what makes a listbox multi-select, and
% what allows an empty selection (needed for "Deselect all").
lst = uicontrol(dlg, 'Style', 'listbox', 'String', itemLabels, 'Units', 'pixels', ...
    'Position', [margin listBottom inner H-52-listBottom], ...
    'Min', 0, 'Max', 2, 'Value', sel, 'FontName', 'Monospaced');

popClass = [];
if hasClasses
    uicontrol(dlg, 'Style', 'text', 'String', 'Group:', 'Units', 'pixels', ...
        'Position', [margin 82 46 18], 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');
    popClass = uicontrol(dlg, 'Style', 'popupmenu', 'String', classNames, ...
        'Units', 'pixels', 'Position', [margin+48 82 142 24]);
    uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Add', 'Units', 'pixels', ...
        'Position', [margin+198 82 74 24], 'Callback', @(s,e) onClassAdd());
    uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Remove', 'Units', 'pixels', ...
        'Position', [margin+278 82 74 24], 'Callback', @(s,e) onClassRemove());
end

uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Select all', 'Units', 'pixels', ...
    'Position', [margin 50 (inner-8)/2 24], 'Callback', @(s,e) set(lst, 'Value', 1:nItems));
uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Deselect all', 'Units', 'pixels', ...
    'Position', [margin+(inner+8)/2 50 (inner-8)/2 24], 'Callback', @(s,e) set(lst, 'Value', []));

uicontrol(dlg, 'Style', 'pushbutton', 'String', confirmStr, 'Units', 'pixels', ...
    'Position', [margin 14 (inner-8)/2 28], 'FontWeight', 'bold', ...
    'Callback', @(s,e) onConfirm());
uicontrol(dlg, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'pixels', ...
    'Position', [margin+(inner+8)/2 14 (inner-8)/2 28], 'Callback', @(s,e) onCancel());

uiwait(dlg);
if isgraphics(dlg), delete(dlg); end


    function onClassAdd()
        % Union, not replace: the point is to seed the selection with a
        % whole category and then hand-tune it.
        set(lst, 'Value', union(get(lst, 'Value'), find(groupMask(:, get(popClass, 'Value')))));
    end

    function onClassRemove()
        set(lst, 'Value', setdiff(get(lst, 'Value'), find(groupMask(:, get(popClass, 'Value')))));
    end

    function onConfirm()
        sel = get(lst, 'Value');
        ok  = true;
        uiresume(dlg);
    end

    function onCancel()
        ok = false;   % sel keeps the incoming preselection
        uiresume(dlg);
    end
end
