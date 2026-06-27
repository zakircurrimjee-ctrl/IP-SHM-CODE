% crack_placement_mse.m
% Modal Strain Energy (MSE) - based optimal crack placement
%
% PURPOSE:
%   Uses the mode shape data saved by Final_EI_STEP.m
%   to compute an MSE-based sensitivity score for every surface node.
%   Nodes with high MSE scores are the best locations for crack placement,
%   because a crack there will produce the largest detectable frequency shift.
%
% THEORY:
%   Modal strain energy of mode i:  MSE_i = 1/2 * phi_i' * K * phi_i (Eq. 3,
%   Shi, Law & Zhang 1998). This is decomposed EXACTLY into a per-node share
%   e_n(i); the per-node score summed over modes ranks each surface location
%   by the eigenvalue (frequency) shift a crack there would produce, since
%   d(lambda_i) = phi_i' * dK * phi_i. Highest score = largest detectable
%   frequency shift. Modes are combined as fractions of their own total energy.
%
% DEPENDENCY:
%      Requires 'ei_workspace.mat' produced by Final_EI_STEP.m
%
% OUTPUT:
%   - Console table of top crack candidate locations (X, Y, Z, MSE score)
%   - Figure 1: MSE heatmap on scattered surface nodes
%   - Figure 2: Per-mode modal strain energy along tube axis
%   - Figure 3: STL model with crack positions overlaid + exported as PNG
%   - crack_candidates.csv       — full ranked node list
%   - crack_positions_on_model.png — publication-ready image of STL + cracks
 
clear; clc; close all;
 
%% =========================
% USER SETTINGS
% ==========================
 
workspaceFile  = 'ei_workspace.mat';   % path to the saved EI workspace
nTopCracks     = 5;                    % how many top crack locations to report
 
% Tube axis: 'auto' detects the longest bounding-box dimension automatically.
% Override with 'x', 'y', or 'z' if auto gives the wrong result.
tubeAxisMode   = 'auto';
 
% ---- Tube geometry (YOUR SPECIFIC TUBE) ----
tubeDiameter_outer = 138.60;  % mm
tubeDiameter_inner = 135.00;  % mm
tubeWallThickness  = (tubeDiameter_outer - tubeDiameter_inner) / 2;  % = 1.80 mm
 
fprintf('Tube geometry:\n');
fprintf('  Outer diameter : %.2f mm\n', tubeDiameter_outer);
fprintf('  Inner diameter : %.2f mm\n', tubeDiameter_inner);
fprintf('  Wall thickness : %.2f mm  <-- very thin, take care with crack depth\n\n', tubeWallThickness);
 
% Minimum axial separation between crack candidates (metres).
% Rule of thumb: ~10–15% of outer diameter keeps candidates physically separated.
% 0.10 * 93.60 mm = ~9.4 mm → rounded to 10 mm.
minAxialSep    = 0.010;   % 10 mm — adjust if your tube is very short
 
exportFigure   = true;   % set false to skip PNG export
 
%% =========================
% 1) LOAD EI WORKSPACE
% ==========================
if ~isfile(workspaceFile)
    error(['Cannot find "%s".\n' ...
       'Run Final_EI_STEP.m first to generate ei_workspace.mat.'], ...
       workspaceFile);
end
 
fprintf('Loading EI workspace from: %s\n', workspaceFile);
load(workspaceFile);
fprintf('Loaded successfully.\n\n');
 
% Validate essential variables
required = {'R','allNodes','freqs_Hz','nModesUse','modeIDs', ...
            'candidateNodeIDs_beforeFilter','candidateXYZ_beforeFilter', ...
            'cadFile','lengthScaleToMeters'};
for k = 1:numel(required)
    if ~exist(required{k}, 'var')
        error('Variable "%s" not found. Re-run the EI script with the save block.', required{k});
    end
end
 
nodeXYZ  = candidateXYZ_beforeFilter;     % [N x 3] surface node coords (metres)
nodeIDs  = candidateNodeIDs_beforeFilter; % [N x 1] global node IDs
N_nodes  = size(nodeXYZ, 1);
N_modes  = nModesUse;
 
fprintf('Surface nodes loaded  : %d\n', N_nodes);
fprintf('Modes used            : %d\n', N_modes);
fprintf('Frequency range       : %.1f Hz  to  %.1f Hz\n\n', ...
        freqs_Hz(1), freqs_Hz(min(N_modes, numel(freqs_Hz))));
 
%% =========================
% 2) DETECT TUBE AXIS
% ==========================
bbox_min  = min(nodeXYZ, [], 1);
bbox_max  = max(nodeXYZ, [], 1);
bbox_span = bbox_max - bbox_min;
 
switch lower(tubeAxisMode)
    case 'auto',  [~, axisIdx] = max(bbox_span);
    case 'x',     axisIdx = 1;
    case 'y',     axisIdx = 2;
    case 'z',     axisIdx = 3;
    otherwise,    error('tubeAxisMode must be ''auto'', ''x'', ''y'', or ''z''.');
end
 
axisNames   = {'X','Y','Z'};
axialCoords = nodeXYZ(:, axisIdx);
 
fprintf('Tube axis detected    : %s  (span = %.4f m)\n\n', ...
        axisNames{axisIdx}, bbox_span(axisIdx));
 
%% =========================
% 2b) FILTER: KEEP ONLY TUBE SURFACE NODES (exclude supports/brackets)
% ==========================
% Nodes on the cylindrical tube wall all sit at the same radial distance
% from the tube axis (= outer radius).  Bracket/support nodes are at
% different radial distances and are excluded by this filter.
%
% Radial distance is computed in the cross-section plane (the two axes
% that are NOT the tube axis).
 
crossIdx = setdiff(1:3, axisIdx);   % e.g. [1 2] when tube axis is Z
 
% Cross-section centroid of ALL candidate nodes
cp = mean(nodeXYZ(:, crossIdx), 1);   % [1 x 2]
 
% Radial distance of every node from the tube axis
r_all = sqrt( (nodeXYZ(:, crossIdx(1)) - cp(1)).^2 + ...
              (nodeXYZ(:, crossIdx(2)) - cp(2)).^2 );   % [N_nodes x 1]
 
% Expected outer radius in metres
tubeOuterRadius_m = (tubeDiameter_outer / 2) * 1e-3;   % 93.60/2 mm → m
 
% Tolerance: nodes within ±radialTol of the outer radius are kept.
% 15% of the outer radius is generous enough to handle mesh irregularities
% while still excluding the flat bracket geometry.
radialTol = tubeOuterRadius_m * 0.15;
 
tubeNodeMask = abs(r_all - tubeOuterRadius_m) <= radialTol;
 
% Apply radial filter
nodeXYZ_all  = nodeXYZ;
nodeIDs_all  = nodeIDs;
nodeXYZ      = nodeXYZ(tubeNodeMask, :);
nodeIDs      = nodeIDs(tubeNodeMask);
axialCoords  = axialCoords(tubeNodeMask);
N_nodes      = size(nodeXYZ, 1);
 
fprintf('Radial filter applied (r = %.2f mm ± %.2f mm):\n', ...
        tubeOuterRadius_m*1e3, radialTol*1e3);
fprintf('  Nodes before filter : %d\n', sum(tubeNodeMask | ~tubeNodeMask));
fprintf('  Nodes after radial  : %d\n', N_nodes);
fprintf('  Removed             : %d\n\n', sum(~tubeNodeMask));
 
%% =========================
% 2c) AXIAL TRIM: EXCLUDE NODES INSIDE THE BRACKET / CLAMP ZONE
% ==========================
% WHY THIS IS NEEDED:
%   Your supports are integrated bracket arms at both ends of the model.
%   The bottom bracket spans roughly Z = 0 to Z = 68 mm.
%   The top bracket spans roughly Z = 274 mm to Z = 290 mm (model top).
%   Nodes in those zones sit on the bracket-tube interface ring — they
%   share the same outer radius as the tube wall (so the radial filter
%   above cannot remove them), but placing a crack there is physically
%   impossible because the bracket clamps and reinforces that zone.
%
% RULE: keep only nodes whose axial position is at least bracketMargin_m
%       away from BOTH ends of the current node set.
%
% Your tube free-wall length is 274 mm. The bracket height visible in
% the model is approximately 68 mm.  A 70 mm margin is therefore safe.
% Reduce this value if it trims too much of the tube.
 
bracketMargin_m = 0.001;  % each end — set to clamped length on your design
 
axialFreeMin = min(axialCoords) + bracketMargin_m;
axialFreeMax = max(axialCoords) - bracketMargin_m;
 
freeWallMask = axialCoords >= axialFreeMin & axialCoords <= axialFreeMax;
 
nodeXYZ     = nodeXYZ(freeWallMask, :);
nodeIDs     = nodeIDs(freeWallMask);
axialCoords = axialCoords(freeWallMask);
N_nodes     = size(nodeXYZ, 1);
 
fprintf('Axial trim applied (bracket margin = %.0f mm each end):\n', bracketMargin_m*1e3);
fprintf('  Free-wall Z range   : %.1f mm  to  %.1f mm\n', ...
        axialFreeMin*1e3, axialFreeMax*1e3);
fprintf('  Nodes in free wall  : %d\n\n', N_nodes);
 
if N_nodes < 10
    error(['Too few tube nodes remain after axial trim (%d). ' ...
           'Reduce bracketMargin_m (currently %.0f mm).'], N_nodes, bracketMargin_m*1e3);
end
 
%% =========================
% 3) MODAL STRAIN ENERGY PER NODE  (true MSE — Shi/Law/Zhang 1998, Eq. 3)
% ==========================
% Replaces the old "curvature of |displacement|" proxy with the actual
% modal strain energy from the FE stiffness matrix.
%
% THEORY
%   Total modal strain energy of mode i:   MSE_i = 1/2 * phi_i' * K * phi_i.
%   This is split EXACTLY into a per-node share:
%       e_n(i) = 1/2 * sum_{c in x,y,z} phi_i[dof(n,c)] * (K*phi_i)[dof(n,c)]
%   so that  sum_n e_n(i) = 1/2 phi_i' K phi_i.
%   A crack at a high-e_n node maximises the eigenvalue shift
%       d(lambda_i) = phi_i' * dK * phi_i,
%   i.e. the frequency change you want to maximise. Modes are combined as
%   fractions of their own total energy so no single mode dominates.
%
% WHY THIS FORM: pde.ModalStructuralResults exposes no strain evaluator, but
% assembleFEMatrices(model) gives the global K. Using K directly needs no
% element shape-function conventions, so it is convention-proof; DOF ordering
% is verified below against the known natural frequency.
 
fprintf('Computing true modal strain energy from the FE stiffness matrix...\n');
 
% --- 3.1 Full mode shapes (all global nodes) ---
if isprop(R.ModeShapes, 'ux')
    ux_all = R.ModeShapes.ux;  uy_all = R.ModeShapes.uy;  uz_all = R.ModeShapes.uz;
elseif isprop(R.ModeShapes, 'x')
    ux_all = R.ModeShapes.x;   uy_all = R.ModeShapes.y;   uz_all = R.ModeShapes.z;
else
    error('Unsupported ModeShapes format. Inspect R.ModeShapes.');
end
nNodesFull = size(ux_all, 1);
UXf = ux_all(:, modeIDs);  UYf = uy_all(:, modeIDs);  UZf = uz_all(:, modeIDs);
 
% --- 3.2 Assemble global stiffness & mass ---
FEM = assembleFEMatrices(model);
K = FEM.K;  M = FEM.M;
if size(K,1) ~= 3*nNodesFull
    error('K is %d but expected 3*nNodes = %d. Unexpected DOF layout.', ...
          size(K,1), 3*nNodesFull);
end
 
% --- 3.3 Detect DOF ordering via the Rayleigh quotient (must equal omega^2) ---
% PDE Toolbox normally interleaves node-major [x1 y1 z1 x2 ...]; we confirm
% against the known frequency so the energy is correct regardless of convention.
buildInterleaved = @(m) reshape([UXf(:,m) UYf(:,m) UZf(:,m)].', 3*nNodesFull, 1);
buildBlock       = @(m) [UXf(:,m); UYf(:,m); UZf(:,m)];
 
wKnown = 2*pi*freqs_Hz(modeIDs(:));   % rad/s
phiI = buildInterleaved(1);  phiB = buildBlock(1);
rqI  = (phiI.'*K*phiI)/(phiI.'*M*phiI);
rqB  = (phiB.'*K*phiB)/(phiB.'*M*phiB);
errI = abs(rqI - wKnown(1)^2)/wKnown(1)^2;
errB = abs(rqB - wKnown(1)^2)/wKnown(1)^2;
if errI <= errB, dofOrder = 'interleaved'; else, dofOrder = 'block'; end
fprintf('  DOF ordering: %s   (Rayleigh error: interleaved=%.2e, block=%.2e)\n', ...
        dofOrder, errI, errB);
if min(errI, errB) > 0.05
    warning(['Rayleigh quotient is >5%% off the known frequency. Energies may ' ...
             'be unreliable — check that mode shapes and K share node numbering.']);
end
 
% --- 3.4 Per-node modal strain energy for each mode ---
e_node   = zeros(nNodesFull, N_modes);
totalMSE = zeros(N_modes, 1);
fprintf('  Per-mode strain energy (verification):\n');
for m = 1:N_modes
    if strcmp(dofOrder,'interleaved'), phi = buildInterleaved(m);
    else,                              phi = buildBlock(m); end
    fK = K*phi;                 % internal force vector
    wk = phi .* fK;             % per-DOF work term (= 2 x energy density)
    if strcmp(dofOrder,'interleaved')
        Wn = reshape(wk, 3, nNodesFull).';                 % [nNodes x 3]
    else
        Wn = [wk(1:nNodesFull), wk(nNodesFull+1:2*nNodesFull), ...
              wk(2*nNodesFull+1:3*nNodesFull)];
    end
    e_node(:,m) = 0.5 * sum(Wn, 2);
    totalMSE(m) = 0.5 * (phi.'*fK);
    relErr = abs(sum(e_node(:,m)) - totalMSE(m)) / max(abs(totalMSE(m)), eps);
    fprintf('    Mode %d (%6.1f Hz):  1/2 phi^T K phi = %.4e  |  expected 1/2 w^2 = %.4e  |  nodal-sum check = %.1e\n', ...
        modeIDs(m), freqs_Hz(modeIDs(m)), totalMSE(m), 0.5*wKnown(m)^2, relErr);
end
 
% --- 3.5 Combine modes as fractions of modal strain energy ---
% Clamp the small negative nodal terms (a discretisation artefact of the nodal
% partition; they appear only in near-zero-energy regions and never affect the
% high-energy ranking). Normalise each mode to a fraction, then sum with equal
% weight (Shi/Law/Zhang recommend combining several modes).
e_pos  = max(e_node, 0);
colSum = sum(e_pos, 1);  colSum(colSum == 0) = 1;
fracN  = e_pos ./ colSum;                 % [nNodesFull x N_modes], each col sums to 1
combinedScore_all = sum(fracN, 2);        % [nNodesFull x 1]
 
% --- 3.6 Map to the placeable surface candidates (after the geometry filters) ---
MSE_score      = combinedScore_all(nodeIDs);     % [N_nodes x 1], aligned to nodeXYZ
MSE_score_norm = MSE_score / max(MSE_score);
fracFiltered   = fracN(nodeIDs, :);              % per-mode, for the axial profile plot
 
fprintf('  Modal strain energy assigned to %d surface candidates.\n\n', N_nodes);
 
%% =========================
% 5) SELECT TOP CRACK CANDIDATES (with axial separation)
% ==========================
[~, rankOrder]  = sort(MSE_score, 'descend');
selectedIdx     = [];
selectedAxial   = [];
 
for k = 1:N_nodes
    candidate = rankOrder(k);
    ax        = axialCoords(candidate);
    if isempty(selectedAxial) || all(abs(selectedAxial - ax) >= minAxialSep)
        selectedIdx(end+1)   = candidate; %#ok<AGROW>
        selectedAxial(end+1) = ax;        %#ok<AGROW>
    end
    if numel(selectedIdx) >= nTopCracks, break; end
end
 
topXYZ    = nodeXYZ(selectedIdx, :);   % metres
topScore  = MSE_score_norm(selectedIdx);
topNodeID = nodeIDs(selectedIdx);
 
%% =========================
% 6) PRINT RESULTS TABLE
% ==========================
fprintf('============================================================\n');
fprintf('  TOP %d CRACK PLACEMENT CANDIDATES (MSE method)\n', nTopCracks);
fprintf('  Tube axis : %s   |   Min axial sep : %.1f mm\n', ...
        axisNames{axisIdx}, minAxialSep*1e3);
fprintf('============================================================\n');
fprintf('  Rank  NodeID   X(mm)    Y(mm)    Z(mm)   MSE score\n');
fprintf('------------------------------------------------------------\n');
for k = 1:numel(selectedIdx)
    fprintf('   %2d   %6d  %7.2f  %7.2f  %7.2f   %.4f\n', ...
        k, topNodeID(k), ...
        topXYZ(k,1)*1e3, topXYZ(k,2)*1e3, topXYZ(k,3)*1e3, ...
        topScore(k));
end
fprintf('============================================================\n');
fprintf('  WALL THICKNESS ADVISORY\n');
fprintf('  Wall = %.2f mm  (OD %.2f mm, ID %.2f mm)\n', ...
        tubeWallThickness, tubeDiameter_outer, tubeDiameter_inner);
fprintf('  Recommended notch depths for SHM dataset diversity:\n');
fprintf('    Shallow  (low severity)  : %.2f mm  (25%% wall)\n', tubeWallThickness * 0.25);
fprintf('    Medium   (mid severity)  : %.2f mm  (50%% wall)\n', tubeWallThickness * 0.50);
fprintf('    Deep     (high severity) : %.2f mm  (75%% wall)\n', tubeWallThickness * 0.75);
fprintf('  Do NOT exceed 75%% depth — risk of fracture during testing.\n');
fprintf('============================================================\n\n');
 
%% =========================
% 7) SAVE CRACK CANDIDATE CSV
% ==========================
crackTable = table( ...
    (1:N_nodes).', nodeIDs(:), ...
    nodeXYZ(:,1)*1e3, nodeXYZ(:,2)*1e3, nodeXYZ(:,3)*1e3, ...
    MSE_score(:), MSE_score_norm(:), ...
    'VariableNames', {'Rank','NodeID','X_mm','Y_mm','Z_mm', ...
                      'MSE_score','MSE_score_normalised'});
crackTable = sortrows(crackTable, 'MSE_score', 'descend');
crackTable.Rank = (1:height(crackTable)).';
writetable(crackTable, 'crack_candidates.csv');
fprintf('Full ranked list saved to: crack_candidates.csv\n\n');
 
%% =========================
% 8) FIGURE 1: MSE HEATMAP ON SCATTERED SURFACE NODES
% ==========================
figure('Name','MSE Sensitivity Map','Color','w');
 
% Show excluded bracket/support nodes in grey for confirmation
excludedXYZ = nodeXYZ_all(~tubeNodeMask, :);
if ~isempty(excludedXYZ)
    scatter3(excludedXYZ(:,1)*1e3, excludedXYZ(:,2)*1e3, excludedXYZ(:,3)*1e3, ...
        10, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.2, ...
        'DisplayName', 'Excluded (supports)');
    hold on;
end
 
% Tube nodes coloured by MSE score
scatter3(nodeXYZ(:,1)*1e3, nodeXYZ(:,2)*1e3, nodeXYZ(:,3)*1e3, ...
    20, MSE_score_norm, 'filled', 'MarkerFaceAlpha', 0.7, ...
    'DisplayName', 'Tube nodes (MSE colour)');
colormap(hot); cb = colorbar;
cb.Label.String = 'Normalised MSE score';
hold on;
 
% Top crack candidates
scatter3(topXYZ(:,1)*1e3, topXYZ(:,2)*1e3, topXYZ(:,3)*1e3, ...
    180, 'cyan', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2, ...
    'DisplayName', 'Top crack candidates');
for k = 1:numel(selectedIdx)
    text(topXYZ(k,1)*1e3+1, topXYZ(k,2)*1e3+1, topXYZ(k,3)*1e3+1, ...
        sprintf(' C%d',k), 'Color','cyan','FontSize',10,'FontWeight','bold', ...
        'BackgroundColor',[0.1 0.1 0.1],'Margin',1);
end
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title({'MSE sensitivity map — tube only (brackets excluded)', ...
       'Hot = high strain energy = best crack location'});
axis equal; grid on; view(3);
legend('Location','best');
 
%% =========================
% 9) FIGURE 2: MODAL STRAIN ENERGY ALONG THE TUBE AXIS
% ==========================
% Axial profile of each mode's strain-energy fraction (binned over the surface
% candidates). Peaks mark the axial bands carrying the most modal strain energy.
figure('Name','Modal strain energy along axis','Color','w');
 
nProf     = 100;
profEdges = linspace(min(axialCoords), max(axialCoords), nProf + 1);
profCtr   = 0.5 * (profEdges(1:end-1) + profEdges(2:end));
profBin   = discretize(axialCoords, profEdges);
profBin(isnan(profBin)) = nProf;
profBin   = min(profBin, nProf);
 
energyProfile = zeros(nProf, N_modes);
for m = 1:N_modes
    for b = 1:nProf
        energyProfile(b,m) = sum(fracFiltered(profBin == b, m));
    end
end
 
cmap = lines(N_modes);
for m = 1:N_modes
    plot(profCtr*1e3, energyProfile(:,m), '-', 'Color', cmap(m,:), ...
        'LineWidth', 1.4, 'DisplayName', ...
        sprintf('Mode %d  (%.1f Hz)', modeIDs(m), freqs_Hz(modeIDs(m))));
    hold on;
end
for k = 1:numel(selectedIdx)
    xline(topXYZ(k,axisIdx)*1e3, '--c', sprintf('C%d',k), ...
        'LabelVerticalAlignment','bottom','LineWidth',1.2);
end
xlabel(sprintf('%s along tube (mm)', axisNames{axisIdx}));
ylabel('Modal strain-energy fraction per axial band');
title('Modal strain energy along tube — peaks = high-sensitivity bands');
legend('Location','best'); grid on;
 
%% =========================
% 10) FIGURE 3: CRACK POSITIONS ON 3D MODEL  (same style as EI sensor plot)
% ==========================
 
% labelOffset scales with tube size so text doesn't overlap the markers.
% ~3% of outer radius works well for your 93.60 mm tube.
labelOffset = (tubeDiameter_outer / 2) * 0.03 * 1e-3;   % in metres
 
fig3 = figure('Name','Crack Positions on 3D Model','Color','w', ...
              'Units','pixels','Position',[100 100 1200 800]);
 
% ---- Base geometry (identical to EI script Section 9) ----
pdegplot(model.Geometry, FaceAlpha=0.10);
hold on;
 
% ---- Excluded bracket/support nodes — shown in grey so you can confirm removal ----
if ~isempty(excludedXYZ)
    scatter3(excludedXYZ(:,1), excludedXYZ(:,2), excludedXYZ(:,3), ...
        10, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15, ...
        'DisplayName', 'Excluded (supports/brackets)');
end
 
% ---- Tube-only surface candidate nodes in green ----
hCandidates = scatter3( ...
    nodeXYZ(:,1), nodeXYZ(:,2), nodeXYZ(:,3), ...
    18, [0.0 0.8 0.0], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25, ...
    'DisplayName', 'Tube surface candidates');
 
% ---- Crack candidates: colour-coded by rank ----
% Colours match the rank order: C1=red, C2=magenta, C3=blue, C4=green, C5=black
% (mirrors the Best/Medium/Worst colouring in the EI plot)
crackMarkerColors = [
    1.00  0.00  0.00;   % C1 — red   (highest MSE)
    0.85  0.00  0.85;   % C2 — magenta
    0.00  0.45  0.90;   % C3 — blue
    0.10  0.70  0.10;   % C4 — green
    0.10  0.10  0.10;   % C5 — black
];
 
nC = numel(selectedIdx);
hCracks = gobjects(nC, 1);
 
for k = 1:nC
    col = crackMarkerColors(mod(k-1, size(crackMarkerColors,1)) + 1, :);
    cx  = topXYZ(k,1);
    cy  = topXYZ(k,2);
    cz  = topXYZ(k,3);
 
    % Filled circle marker — same size as sensor markers in EI script
    hCracks(k) = scatter3(cx, cy, cz, 120, col, 'filled', ...
        'DisplayName', sprintf('C%d  (%s=%.1f mm, MSE=%.3f)', ...
        k, axisNames{axisIdx}, topXYZ(k,axisIdx)*1e3, topScore(k)));
 
    % Direction arrow pointing radially outward from tube axis
    % (mirrors the quiver3 sensor direction arrows in the EI script)
    tubeCenter = mean(nodeXYZ(:, setdiff(1:3, axisIdx)), 1);  % centroid in cross-section
    radialVec  = [cx cy cz];
    radialVec(axisIdx) = tubeCenter(axisIdx - (axisIdx>1));
    radialVec  = [cx cy cz] - [tubeCenter(1) tubeCenter(1) tubeCenter(1)];
 
    % Robust radial direction: point from tube axis to crack node
    crossIdx   = setdiff(1:3, axisIdx);   % the two non-axial dimensions
    cp         = mean(nodeXYZ(:, crossIdx), 1);   % tube cross-section centroid
    rv         = zeros(1,3);
    nodeVec      = [cx cy cz];
    rv(crossIdx) = nodeVec(crossIdx) - cp;
    rv_norm    = rv / (norm(rv) + eps);
 
    arrowScale = (tubeDiameter_outer/2) * 0.25 * 1e-3;   % 25% of outer radius
    quiver3(cx, cy, cz, ...
            rv_norm(1)*arrowScale, rv_norm(2)*arrowScale, rv_norm(3)*arrowScale, ...
            0, 'Color', col, 'LineWidth', 1.8, 'HandleVisibility', 'off');
 
    % Text label — same style as B1/M1/W1 labels in EI script
    text(cx + labelOffset, cy + labelOffset, cz + labelOffset, ...
        sprintf('C%d', k), ...
        'Color',           col, ...
        'FontSize',        10, ...
        'FontWeight',      'bold', ...
        'BackgroundColor', 'w', ...
        'Margin',          1);
end
 
% ---- Axes, title, legend ----
legend([hCandidates; hCracks], 'Location', 'best');
 
title(sprintf('MSE-based crack placement — top %d candidates (tube axis: %s)', ...
      nC, axisNames{axisIdx}));
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal; grid on; view(3);
 
% ---- Export high-resolution PNG ----
if exportFigure
    exportFile = 'crack_positions_on_model.png';
    exportgraphics(fig3, exportFile, 'Resolution', 300);
    fprintf('Figure exported to: %s\n\n', exportFile);
end
 
fprintf('All done.\n');
fprintf('  crack_candidates.csv         — full ranked node list\n');
if exportFigure
    fprintf('  crack_positions_on_model.png — 3D model with crack markers\n');
end