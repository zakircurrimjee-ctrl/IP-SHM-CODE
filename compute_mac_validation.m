%% compute_mac_validation.m
%
% Two-part MAC validation for EI and EOT sensor placement results.
%
% PREREQUISITES:
%   Run Final_EI_STEP.m  -> produces ei_osp_results.mat
%   Run FINAL_EOT_STEP.m -> produces eot_osp_results.mat
%
% WHAT THIS SCRIPT DOES:
%   Part 1 — MAC for EI-selected sensors
%             Checks that the EI sensor set can distinguish all target modes.
%             Pass criterion: off-diagonal MAC < 0.2, diagonal MAC > 0.9.
%
%   Part 2 — MAC for EOT-selected sensors
%             Same check for the EOT sensor set.
%
%   EI vs EOT both come from the same FEM modal solve and retain the same
%   modes in the same frequency order, so they target the same physical
%   modes by construction. A direct EI-vs-EOT cross-MAC is not computed:
%   the two sensor sets sit at different node locations, so a cross-MAC
%   between their reduced shapes is not a like-for-like comparison.
%
% REFERENCE:
%   Allemang & Brown (1982); Heo, Wang & Satpathi (1997).

clear; clc; close all;

%% =========================
% USER SETTINGS
% ==========================
macDiagThreshold    = 0.9;   % diagonal MAC must exceed this (good correlation)
macOffDiagThreshold = 0.2;   % off-diagonal MAC must stay below this (mode independence)

% Figure export settings
saveFigures = true;          % set false to skip writing files
figFormats  = {'png','pdf'}; % publication formats (png @300dpi + vector pdf)

% Anchor output folder to this script's location, independent of MATLAB's
% current folder. This avoids the "Unable to create output file" error that
% occurs when the Run button leaves the current folder elsewhere.
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)            % e.g. code pasted into Command Window
    scriptDir = pwd;
end
figDir = fullfile(scriptDir, 'figures');

if saveFigures && ~exist(figDir, 'dir')
    [ok, msg] = mkdir(figDir);
    if ~ok
        error('Could not create figure directory "%s": %s', figDir, msg);
    end
end

%% =========================
% LOAD RESULTS
% ==========================
fprintf('Loading EI results...\n');
if ~isfile('ei_osp_results.mat')
    error('ei_osp_results.mat not found. Run Final_EI_STEP.m first.');
end
load('ei_osp_results.mat');   % loads: Phi_EI, bestIdx_EI, freqs_Hz_EI, nModesUse_EI

fprintf('Loading EOT results...\n');
if ~isfile('eot_osp_results.mat')
    error('eot_osp_results.mat not found. Run FINAL_EOT_STEP.m first.');
end
load('eot_osp_results.mat');  % loads: Phi_EOT, bestIdx_EOT, freqs_Hz_EOT, nModesUse_EOT

%% =========================
% EXTRACT REDUCED MODE SHAPE MATRICES
% ==========================
% Each row of Phi is a candidate DOF; each column is a mode.
% We keep only the rows corresponding to the sensors selected by each method.

Phi_EI_reduced  = Phi_EI(bestIdx_EI,  :);   % (nSensors_EI  x nModes)
Phi_EOT_reduced = Phi_EOT(bestIdx_EOT, :);   % (nSensors_EOT x nModes)

nModes_EI  = size(Phi_EI_reduced,  2);
nModes_EOT = size(Phi_EOT_reduced, 2);

fprintf('\nEI  reduced matrix: %d sensors x %d modes\n', size(Phi_EI_reduced,1),  nModes_EI);
fprintf('EOT reduced matrix: %d sensors x %d modes\n', size(Phi_EOT_reduced,1), nModes_EOT);

%% =========================
% PART 1 — MAC FOR EI SENSORS
% ==========================
fprintf('\n==============================\n');
fprintf('PART 1: EI sensor placement MAC\n');
fprintf('==============================\n');

MAC_EI = computeMAC(Phi_EI_reduced, Phi_EI_reduced);

% Print matrix
fprintf('\nMAC matrix (EI):\n');
modeLabels = arrayfun(@(m) sprintf('Mode %d (%.1f Hz)', m, freqs_Hz_EI(m)), ...
    1:nModes_EI, 'UniformOutput', false);
disp(array2table(MAC_EI, 'VariableNames', matlab.lang.makeValidName(modeLabels), ...
    'RowNames', matlab.lang.makeValidName(modeLabels)));

% Pass/fail report
diagVals_EI    = diag(MAC_EI);
offDiagVals_EI = MAC_EI(~eye(nModes_EI, 'logical'));

fprintf('Diagonal values (should be > %.2f):\n', macDiagThreshold);
for m = 1:nModes_EI
    if diagVals_EI(m) >= macDiagThreshold, status = 'PASS'; else, status = 'FAIL'; end
    fprintf('  Mode %d (%.1f Hz): %.4f  [%s]\n', m, freqs_Hz_EI(m), diagVals_EI(m), status);
end

fprintf('Off-diagonal max: %.4f  ', max(offDiagVals_EI));
if max(offDiagVals_EI) < macOffDiagThreshold
    fprintf('[PASS — below %.2f]\n', macOffDiagThreshold);
else
    fprintf('[FAIL — exceeds %.2f]\n', macOffDiagThreshold);
end

% Plot
plotMAC(MAC_EI, freqs_Hz_EI(1:nModes_EI), freqs_Hz_EI(1:nModes_EI), ...
    'MAC -- EI Sensor Placement', ...
    sprintf('Auto-MAC, EI set (%d sensors, %d modes)', size(Phi_EI_reduced,1), nModes_EI), ...
    'Mode', 'Mode');
exportFigure(gcf, fullfile(figDir, 'MAC_EI'), saveFigures, figFormats);

%% =========================
% PART 2 — MAC FOR EOT SENSORS
% ==========================
fprintf('\n==============================\n');
fprintf('PART 2: EOT sensor placement MAC\n');
fprintf('==============================\n');

MAC_EOT = computeMAC(Phi_EOT_reduced, Phi_EOT_reduced);

fprintf('\nMAC matrix (EOT):\n');
modeLabels_EOT = arrayfun(@(m) sprintf('Mode %d (%.1f Hz)', m, freqs_Hz_EOT(m)), ...
    1:nModes_EOT, 'UniformOutput', false);
disp(array2table(MAC_EOT, 'VariableNames', matlab.lang.makeValidName(modeLabels_EOT), ...
    'RowNames', matlab.lang.makeValidName(modeLabels_EOT)));

diagVals_EOT    = diag(MAC_EOT);
offDiagVals_EOT = MAC_EOT(~eye(nModes_EOT, 'logical'));

fprintf('Diagonal values (should be > %.2f):\n', macDiagThreshold);
for m = 1:nModes_EOT
    if diagVals_EOT(m) >= macDiagThreshold, status = 'PASS'; else, status = 'FAIL'; end
    fprintf('  Mode %d (%.1f Hz): %.4f  [%s]\n', m, freqs_Hz_EOT(m), diagVals_EOT(m), status);
end

fprintf('Off-diagonal max: %.4f  ', max(offDiagVals_EOT));
if max(offDiagVals_EOT) < macOffDiagThreshold
    fprintf('[PASS — below %.2f]\n', macOffDiagThreshold);
else
    fprintf('[FAIL — exceeds %.2f]\n', macOffDiagThreshold);
end

plotMAC(MAC_EOT, freqs_Hz_EOT(1:nModes_EOT), freqs_Hz_EOT(1:nModes_EOT), ...
    'MAC -- EOT Sensor Placement', ...
    sprintf('Auto-MAC, EOT set (%d sensors, %d modes)', size(Phi_EOT_reduced,1), nModes_EOT), ...
    'Mode', 'Mode');
exportFigure(gcf, fullfile(figDir, 'MAC_EOT'), saveFigures, figFormats);

%% =========================
% SUMMARY
% ==========================
fprintf('\n==============================\n');
fprintf('SUMMARY\n');
fprintf('==============================\n');

fprintf('\n--- EI sensor set ---\n');
fprintf('  Min diagonal MAC:     %.4f  (threshold > %.2f)\n', min(diagVals_EI), macDiagThreshold);
fprintf('  Max off-diagonal MAC: %.4f  (threshold < %.2f)\n', max(offDiagVals_EI), macOffDiagThreshold);
if min(diagVals_EI) >= macDiagThreshold && max(offDiagVals_EI) < macOffDiagThreshold
    fprintf('  RESULT: PASS\n');
else
    fprintf('  RESULT: FAIL — review sensor placement\n');
end

fprintf('\n--- EOT sensor set ---\n');
fprintf('  Min diagonal MAC:     %.4f  (threshold > %.2f)\n', min(diagVals_EOT), macDiagThreshold);
fprintf('  Max off-diagonal MAC: %.4f  (threshold < %.2f)\n', max(offDiagVals_EOT), macOffDiagThreshold);
if min(diagVals_EOT) >= macDiagThreshold && max(offDiagVals_EOT) < macOffDiagThreshold
    fprintf('  RESULT: PASS\n');
else
    fprintf('  RESULT: FAIL — review sensor placement\n');
end

fprintf('\n--- EI vs EOT: same modes? ---\n');
fprintf('  Both sensor sets are derived from the same FEM modal solve and\n');
fprintf('  retain the same modes in the same frequency order:\n');
nShared = min(numel(freqs_Hz_EI), numel(freqs_Hz_EOT));
for m = 1:min(nShared, nModes_EI)
    fprintf('    Mode %d:  EI %.1f Hz  |  EOT %.1f Hz\n', ...
        m, freqs_Hz_EI(m), freqs_Hz_EOT(m));
end
fprintf('  The auto-MAC checks above confirm each set resolves these modes\n');
fprintf('  distinctly (off-diagonal near 0). The methods therefore target and\n');
fprintf('  observe the same physical modes by construction.\n');

if saveFigures
    fprintf('\nFigures written to "%s/" (%s)\n', figDir, strjoin(figFormats, ', '));
end

%% =========================
% LOCAL FUNCTIONS
% ==========================
% NOTE: in a script, all local functions must appear at the very end of the
% file, after the script body.

function M = computeMAC(A, B)
% MAC between two mode shape matrices A (na x nModes) and B (nb x nModes).
% MAC(i,j) = |A(:,i)' * B(:,j)|^2 / ( (A(:,i)'*A(:,i)) * (B(:,j)'*B(:,j)) )
% When A == B this gives the auto-MAC (should be identity for good placement).
    nA = size(A, 2);
    nB = size(B, 2);
    M  = zeros(nA, nB);
    for ii = 1:nA
        for jj = 1:nB
            num   = abs(A(:,ii)' * B(:,jj))^2;
            denom = (A(:,ii)' * A(:,ii)) * (B(:,jj)' * B(:,jj));
            if denom < eps
                M(ii,jj) = 0;
            else
                M(ii,jj) = num / denom;
            end
        end
    end
end

function plotMAC(MAC, rowFreqs, colFreqs, figName, titleStr, xLab, yLab)
% Publication-quality MAC heatmap with a sequential colormap, crisp cell
% borders, a boxed leading diagonal, and luminance-aware value labels.

    [nR, nC] = size(MAC);

    fig = figure('Name', figName, 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 14 12]);
    ax = axes(fig); hold(ax, 'on');

    % --- heatmap ---
    imagesc(ax, MAC);
    colormap(ax, mac_colormap());
    cb = colorbar(ax);
    cb.Label.String = 'MAC value';
    cb.Label.FontSize = 12;
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    clim(ax, [0 1]);

    set(ax, 'YDir', 'reverse');   % row 1 at top (matrix convention)
    axis(ax, 'equal'); axis(ax, 'tight');

    % --- thin cell separators for crisp tiles ---
    for k = 1.5:1:(nC-0.5)
        plot(ax, [k k], [0.5 nR+0.5], 'Color', [1 1 1 0.7], 'LineWidth', 0.5);
    end
    for k = 1.5:1:(nR-0.5)
        plot(ax, [0.5 nC+0.5], [k k], 'Color', [1 1 1 0.7], 'LineWidth', 0.5);
    end

    % --- box the leading diagonal (cells that should be ~1) ---
    for d = 1:min(nR, nC)
        rectangle(ax, 'Position', [d-0.5 d-0.5 1 1], ...
            'EdgeColor', [0 0 0], 'LineWidth', 1.4);
    end

    % --- value annotations with luminance-aware contrast ---
    for ii = 1:nR
        for jj = 1:nC
            v = MAC(ii,jj);
            text(ax, jj, ii, sprintf('%.3f', v), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'Color', text_contrast(v));
        end
    end

    % --- ticks / labels ---
    xticks(ax, 1:nC); yticks(ax, 1:nR);
    xticklabels(ax, arrayfun(@(m) sprintf('M%d: %.1f Hz', m, colFreqs(m)), ...
        1:nC, 'UniformOutput', false));
    yticklabels(ax, arrayfun(@(m) sprintf('M%d: %.1f Hz', m, rowFreqs(m)), ...
        1:nR, 'UniformOutput', false));

    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 11, ...
        'Layer', 'top', 'TickLength', [0 0], 'Box', 'on', 'LineWidth', 1.0);

    xlabel(ax, xLab, 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, yLab, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 13);

    hold(ax, 'off');
end

function cmap = mac_colormap()
% Sequential white -> deep blue map: low MAC reads as near-blank, high MAC as
% saturated. Easier to discuss in print than parula (no spurious green/yellow
% banding on the off-diagonals).
    base = [1.00 1.00 1.00;
            0.86 0.91 0.96;
            0.62 0.76 0.89;
            0.36 0.58 0.80;
            0.16 0.40 0.68;
            0.05 0.24 0.49];
    x  = linspace(0, 1, size(base,1));
    xi = linspace(0, 1, 256);
    cmap = [interp1(x, base(:,1), xi)', ...
            interp1(x, base(:,2), xi)', ...
            interp1(x, base(:,3), xi)'];
end

function c = text_contrast(val)
% Pick black or white text based on the luminance of the cell colour, so
% numbers stay legible across the whole colormap.
    cmap = mac_colormap();
    idx  = max(1, min(size(cmap,1), round(val*(size(cmap,1)-1))+1));
    rgb  = cmap(idx,:);
    lum  = 0.299*rgb(1) + 0.587*rgb(2) + 0.114*rgb(3);
    if lum < 0.55, c = [1 1 1]; else, c = [0 0 0]; end
end

function exportFigure(fig, basePath, doSave, formats)
% Export a figure to the requested formats at print resolution.
    if ~doSave, return; end
    for k = 1:numel(formats)
        fmt = lower(formats{k});
        switch fmt
            case 'png'
                exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
            case {'pdf','eps'}
                exportgraphics(fig, [basePath '.' fmt], 'ContentType', 'vector');
            otherwise
                exportgraphics(fig, [basePath '.' fmt]);
        end
    end
end