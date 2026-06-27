%% Map Hammer Strike Positions onto the STEP Geometry and the Physical Tube
% --------------------------------------------------------------------------
% PURPOSE
%   Take the ADPR-ranked strike points and (a) overlay them on the actual CAD
%   surface render, and (b) convert each to bench-usable coordinates so you
%   can mark the real specimen with a ruler and protractor:
%       - axial distance measured from the tube end (Z = 0 collar face), and
%       - circumferential angle around the tube.
%
% GEOMETRY ASSUMPTION (STEP file)
%   Tube axis runs along GLOBAL Z (Z extent approx -0.299 .. 0 m, outer
%   radius approx 0.0693 m from the Z axis). The angle is therefore measured
%   in the XY plane.
%
% ANGULAR CONVENTION
%   theta = atan2(Y, X) measured from the +X axis, right-handed about +Z,
%   reported in degrees in [0, 360). If you prefer +Y as zero, see the
%   angleZeroAxis option below.
% --------------------------------------------------------------------------

clear; clc; close all;

%% =========================
% 0) LOAD + RECOMPUTE ADPR RANKING
% ==========================
% This script is self-contained: it reloads the workspace and recomputes the
% ADPR ranking using the same surface-normal measurement basis as the OSP scripts
% on its own. Keep the two metric definitions in sync if you edit one.

wsFile = "ei_workspace.mat";
if ~isfile(wsFile)
   error('%s not found. Run Final_EI_STEP.m first.', wsFile);
end
S = load(wsFile);

R                = S.R;
allNodes         = S.allNodes;
candidateNodeIDs = S.candidateNodeIDs_beforeFilter;
candidateXYZ     = S.candidateXYZ_beforeFilter;   % meters
mesh             = S.mesh;                         % needed for surface-normal sensor dirs
modeIDs          = S.modeIDs;
freqs_Hz         = S.freqs_Hz;
hasModel         = isfield(S, 'model');
if hasModel, model = S.model; end

nModes  = numel(modeIDs);

%% =========================
% 1) OPTIONS
% ==========================
% Tube axis for STEP geometry = global Z.
axisPoint1   = [0 0 0];
axisPoint2   = [0 0 1];
axisVec      = (axisPoint2 - axisPoint1); axisVec = axisVec / norm(axisVec);

nReport      = 10;          % strike points to map
angleZeroAxis = "X";        % "X" -> 0 deg at +X ; "Y" -> 0 deg at +Y
axialFromEnd  = "Zmax";     % "Zmax" -> distance measured from Z = max(Z) end
                            % "Zmin" -> distance measured from Z = min(Z) end

%% =========================
% 2) ADPR RANKING (same metric as drive_point_adpr.m)
% ==========================
% Sensor measurement direction = local outward surface normal, matching the
% new surfaceNormal-only EI/EOT scripts (Final_EI_STEP.m / FINAL_EOT_STEP.m).
% This replaces the previous radialDirectionsFromAxis projection so the ADPR
% ranking uses the same measurement basis the OSP scripts actually used.
sensorDirs = surfaceNormalsAtNodesFromTetMesh(mesh, allNodes, candidateNodeIDs);
Phi    = buildProjectedModalMatrix(R, candidateNodeIDs, sensorDirs, modeIDs);
absPhi = abs(Phi);
colMax = max(absPhi, [], 1); colMax(colMax == 0) = 1;
absPhiN = absPhi ./ colMax;
ADPR   = exp( mean( log(absPhiN + 1e-12), 2 ) );

[~, order] = sort(ADPR, 'descend');
nReport = min(nReport, numel(order));
topIdx  = order(1:nReport);

%% =========================
% 3) SELF-CHECK THE AXIS ASSUMPTION
% ==========================
% Radial distance of every candidate node from the assumed axis. On the outer
% wall these should all be ~constant (the outer radius). If they are wildly
% varying, the axis is wrong and all angles/axial distances are meaningless.
rel        = candidateXYZ - axisPoint1;
axialProj  = (rel * axisVec(:));                 % scalar projection onto axis
radialVec  = rel - axialProj * axisVec;          % component perpendicular to axis
radii      = vecnorm(radialVec, 2, 2);

rMean = mean(radii); rStd = std(radii);
fprintf('--- Axis self-check (assumed axis = Z) ---\n');
fprintf('  radial distance from axis: mean = %.2f mm, std = %.2f mm\n', ...
    rMean*1000, rStd*1000);
if rStd / max(rMean, eps) > 0.05
    warning(['Radii vary by >5%% of the mean. The assumed Z axis may be ', ...
             'WRONG, or the candidate face spans collars of different radius. ', ...
             'Angles and axial distances below may be unreliable. Inspect the ', ...
             'overlay figure before trusting them.']);
else
    fprintf('  OK: radii ~constant, Z-axis assumption is consistent.\n');
end
fprintf('\n');

%% =========================
% 4) CONVERT TO BENCH COORDINATES (axial + angle)
% ==========================
% Axial coordinate along the tube.
Zc   = candidateXYZ(:,3);
switch upper(axialFromEnd)
    case "ZMAX", axialDist = max(Zc) - Zc;   % distance from the Z = max end
    case "ZMIN", axialDist = Zc - min(Zc);   % distance from the Z = min end
    otherwise,   error('axialFromEnd must be "Zmax" or "Zmin".');
end

% Angle in the XY plane.
switch upper(angleZeroAxis)
    case "X", theta = atan2d(candidateXYZ(:,2), candidateXYZ(:,1));
    case "Y", theta = atan2d(-candidateXYZ(:,1), candidateXYZ(:,2));
    otherwise, error('angleZeroAxis must be "X" or "Y".');
end
theta = mod(theta, 360);   % wrap to [0, 360)

%% =========================
% 5) REPORT
% ==========================
fprintf('=== Strike points mapped to specimen (top %d) ===\n', nReport);
fprintf('Axial distance measured from the Z = %s end of the tube.\n', upper(axialFromEnd));
fprintf('Angle measured from +%s axis, right-handed about +Z.\n\n', upper(angleZeroAxis));
fprintf('%4s  %8s  %9s  %9s  %9s  %10s  %8s  %8s\n', ...
    'rank','nodeID','X(mm)','Y(mm)','Z(mm)','axial(mm)','angle','ADPR');
for k = 1:nReport
    i  = topIdx(k);
    id = candidateNodeIDs(i);
    xyz = candidateXYZ(i,:)*1000;
    fprintf('%4d  %8d  %9.2f  %9.2f  %9.2f  %10.2f  %7.1f°  %8.4f\n', ...
        k, id, xyz(1), xyz(2), xyz(3), axialDist(i)*1000, theta(i), ADPR(i));
end

best = topIdx(1);
fprintf('\nBEST STRIKE POINT on the physical tube:\n');
fprintf('  %.1f mm from the Z=%s end, at %.1f° around the tube.\n', ...
    axialDist(best)*1000, upper(axialFromEnd), theta(best));
fprintf('  (global XYZ: [%.2f, %.2f, %.2f] mm)\n', candidateXYZ(best,:)*1000);

%% =========================
% 5b) PER-MODE vs COMBINED BEST STRIKE (modes 1 and 2 side by side)
% ==========================
% Three "best" points:
%   - mode 1 alone : node with the largest |Phi(:,1)| (most excitable for mode 1)
%   - mode 2 alone : node with the largest |Phi(:,2)|
%   - combined     : the ADPR winner (best compromise that excites BOTH)
% Per-mode ranking uses |Phi(:,m)|, the single-mode analogue of ADPR. The three
% points generally differ: each single mode peaks where its own antinode sits,
% while the combined point is pulled toward a location good for both at once.

% Index (into the candidate list) of the best node for each individual mode.
[~, bestMode1] = max(absPhi(:,1));
[~, bestMode2] = max(absPhi(:,2));
bestCombined   = best;   % ADPR winner from Section 4

pickRows = [bestMode1, bestMode2, bestCombined];
pickName = ["Mode 1 only"; "Mode 2 only"; "Combined (ADPR)"];

fprintf('\n=== Best strike: per-mode vs combined ===\n');
fprintf('Mode 1 = %.2f Hz, Mode 2 = %.2f Hz.\n', ...
    freqs_Hz(modeIDs(1)), freqs_Hz(modeIDs(2)));
fprintf('"|phi1|" and "|phi2|" are the normalized excitabilities (0..1) of each\n');
fprintf('mode AT that point: 1 = that mode''s antinode, ~0 = a node for that mode.\n\n');
fprintf('%-16s  %8s  %10s  %8s  %8s  %8s  %8s\n', ...
    'target','nodeID','axial(mm)','angle','|phi1|','|phi2|','ADPR');
for r = 1:numel(pickRows)
    i  = pickRows(r);
    id = candidateNodeIDs(i);
    fprintf('%-16s  %8d  %10.2f  %7.1f°  %8.3f  %8.3f  %8.4f\n', ...
        pickName(r), id, axialDist(i)*1000, theta(i), ...
        absPhiN(i,1), absPhiN(i,2), ADPR(i));
end

fprintf('\nReading this table:\n');
fprintf('  - For exciting ONE mode in isolation, use its own row (highest |phi| for that mode).\n');
fprintf('  - For one hit that captures BOTH modes, use the Combined row.\n');
fprintf('  - If the Combined row already has high |phi1| AND |phi2|, a single strike\n');
fprintf('    there is enough; no need for separate per-mode hits.\n');

%% =========================
% 6) OVERLAY ON CAD SURFACE  (mirrors the working EI script figure)
% ==========================
% Define best / medium / worst strike tiers from the ADPR ranking, the same
% way the EI script shows best/medium/worst SENSORS. Here "best" = excites all
% modes well (high ADPR), "worst" = lands near a nodal line (low ADPR).
nTier = min(3, numel(order));               % markers per tier
bestStrikeIdx   = order(1:nTier);                        % top ADPR
worstStrikeIdx  = order(end-nTier+1:end);               % bottom ADPR
midStart        = max(1, round(numel(order)/2) - floor(nTier/2));
mediumStrikeIdx = order(midStart:midStart+nTier-1);     % middle ADPR

bestStrikeXYZ   = candidateXYZ(bestStrikeIdx,:);
mediumStrikeXYZ = candidateXYZ(mediumStrikeIdx,:);
worstStrikeXYZ  = candidateXYZ(worstStrikeIdx,:);

figure('Color','w','Name','Strike points on CAD surface');

% IMPORTANT: draw the geometry FIRST, then hold on. Calling hold on before
% pdegplot makes the geometry render fight the auto-ranging and collapses the
% view to a tiny cluster. This order matches ei_osp_fast_uprgradedSTEP.m.
if hasModel
    pdegplot(model.Geometry, FaceAlpha=0.10);
else
    scatter3(candidateXYZ(:,1), candidateXYZ(:,2), candidateXYZ(:,3), ...
        4, [0.8 0.8 0.8], 'filled');
end
hold on;

% Faint candidate cloud, colored by ADPR for context.
hCand = scatter3(candidateXYZ(:,1), candidateXYZ(:,2), candidateXYZ(:,3), ...
    18, ADPR, 'filled', 'MarkerFaceAlpha', 0.35);

% Tier markers: best = red, medium = magenta, worst = black (same palette
% as the EI sensor figure so the two are visually consistent).
hBest = scatter3(bestStrikeXYZ(:,1), bestStrikeXYZ(:,2), bestStrikeXYZ(:,3), ...
    140, 'r', 'filled', 'MarkerEdgeColor','k');
hMedium = scatter3(mediumStrikeXYZ(:,1), mediumStrikeXYZ(:,2), mediumStrikeXYZ(:,3), ...
    140, 'm', 'filled', 'MarkerEdgeColor','k');
hWorst = scatter3(worstStrikeXYZ(:,1), worstStrikeXYZ(:,2), worstStrikeXYZ(:,3), ...
    140, 'k', 'filled');

% Strike-direction arrows (surface-normal outward), same style as the sensor quivers.
bestStrikeDirs = sensorDirs(bestStrikeIdx,:);
quiver3(bestStrikeXYZ(:,1), bestStrikeXYZ(:,2), bestStrikeXYZ(:,3), ...
        bestStrikeDirs(:,1), bestStrikeDirs(:,2), bestStrikeDirs(:,3), ...
        0.03, 'r', 'LineWidth', 1.5, 'HandleVisibility','off');

% Labels with white background, matching the EI figure (B/M/W tiers).
labelOffset = 0.003;   % meters
for i = 1:size(bestStrikeXYZ,1)
    text(bestStrikeXYZ(i,1)+labelOffset, bestStrikeXYZ(i,2)+labelOffset, bestStrikeXYZ(i,3)+labelOffset, ...
        sprintf('B%d', i), 'Color','r', 'FontSize',10, 'FontWeight','bold', ...
        'BackgroundColor','w', 'Margin',1);
end
for i = 1:size(mediumStrikeXYZ,1)
    text(mediumStrikeXYZ(i,1)+labelOffset, mediumStrikeXYZ(i,2)+labelOffset, mediumStrikeXYZ(i,3)+labelOffset, ...
        sprintf('M%d', i), 'Color','m', 'FontSize',10, 'FontWeight','bold', ...
        'BackgroundColor','w', 'Margin',1);
end
for i = 1:size(worstStrikeXYZ,1)
    text(worstStrikeXYZ(i,1)+labelOffset, worstStrikeXYZ(i,2)+labelOffset, worstStrikeXYZ(i,3)+labelOffset, ...
        sprintf('W%d', i), 'Color','k', 'FontSize',10, 'FontWeight','bold', ...
        'BackgroundColor','w', 'Margin',1);
end

colormap(parula); cb = colorbar; cb.Label.String = 'ADPR (higher = better)';
legend([hCand, hBest, hMedium, hWorst], ...
    {'Candidate nodes (ADPR)', 'Best strike', 'Medium strike', 'Worst strike'}, ...
    'Location', 'best');
title(sprintf('Best / Medium / Worst hammer strike positions (%d modes)', nModes));
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal; grid on; view(3);
hold off;

%% =========================
% 7) PER-MODE vs COMBINED STRIKE POINTS (3-point figure)
% ==========================
% Three points from Section 5b: best for mode 1 alone, best for mode 2 alone,
% and the combined ADPR point. Same draw order as above (geometry first, then
% hold on) so the tube renders solid and the view frames correctly.
xyz1 = candidateXYZ(bestMode1,:);    % best for mode 1
xyz2 = candidateXYZ(bestMode2,:);    % best for mode 2
xyzC = candidateXYZ(bestCombined,:); % combined (ADPR)

figure('Color','w','Name','Per-mode vs combined strike points');

if hasModel
    pdegplot(model.Geometry, FaceAlpha=0.10);
else
    scatter3(candidateXYZ(:,1), candidateXYZ(:,2), candidateXYZ(:,3), ...
        4, [0.85 0.85 0.85], 'filled');
end
hold on;

% Faint candidate cloud for context (single neutral color, not ADPR, so the
% three points stand out).
scatter3(candidateXYZ(:,1), candidateXYZ(:,2), candidateXYZ(:,3), ...
    14, [0.75 0.75 0.85], 'filled', 'MarkerFaceAlpha', 0.30, ...
    'HandleVisibility','off');

% The three points shown as ARROWS pointing inward to each spot. Each arrow
% starts a short distance radially OUTSIDE the wall and points in to the node,
% so the arrow tip marks the exact strike location ("hit here"). Colors:
% mode 1 = blue, mode 2 = green, combined = red.
arrowLen = 0.04;   % meters; visible length of each arrow
d1 = sensorDirs(bestMode1,:);    % outward surface-normal unit vectors
d2 = sensorDirs(bestMode2,:);
dC = sensorDirs(bestCombined,:);

% Tail = point + outward*arrowLen ; direction = inward (-outward), so tip lands
% on the node.
quiver3(xyz1(1)+d1(1)*arrowLen, xyz1(2)+d1(2)*arrowLen, xyz1(3)+d1(3)*arrowLen, ...
        -d1(1), -d1(2), -d1(3), arrowLen, 'b', 'LineWidth', 2.5, ...
        'MaxHeadSize', 1.0, 'HandleVisibility','off');
quiver3(xyz2(1)+d2(1)*arrowLen, xyz2(2)+d2(2)*arrowLen, xyz2(3)+d2(3)*arrowLen, ...
        -d2(1), -d2(2), -d2(3), arrowLen, 'g', 'LineWidth', 2.5, ...
        'MaxHeadSize', 1.0, 'HandleVisibility','off');
quiver3(xyzC(1)+dC(1)*arrowLen, xyzC(2)+dC(2)*arrowLen, xyzC(3)+dC(3)*arrowLen, ...
        -dC(1), -dC(2), -dC(3), arrowLen, 'r', 'LineWidth', 2.5, ...
        'MaxHeadSize', 1.0, 'HandleVisibility','off');

% Small dots exactly at each tip so the landing point is unambiguous, and to
% carry the legend entries (quiver3 does not legend cleanly).
h1 = scatter3(xyz1(1), xyz1(2), xyz1(3), 30, 'b', 'filled');
h2 = scatter3(xyz2(1), xyz2(2), xyz2(3), 30, 'g', 'filled');
hC = scatter3(xyzC(1), xyzC(2), xyzC(3), 40, 'r', 'filled');

% Labels.
lo = 0.004;
text(xyz1(1)+lo, xyz1(2)+lo, xyz1(3)+lo, 'Mode 1', 'Color','b', ...
    'FontSize',10, 'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
text(xyz2(1)+lo, xyz2(2)+lo, xyz2(3)+lo, 'Mode 2', 'Color',[0 0.5 0], ...
    'FontSize',10, 'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
text(xyzC(1)+lo, xyzC(2)+lo, xyzC(3)+lo, 'Both', 'Color','r', ...
    'FontSize',10, 'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);

legend([h1, h2, hC], ...
    {sprintf('Best for mode 1 (%.1f Hz)', freqs_Hz(modeIDs(1))), ...
     sprintf('Best for mode 2 (%.1f Hz)', freqs_Hz(modeIDs(2))), ...
     'Best for both (ADPR)'}, ...
    'Location', 'best');
title('Best strike point: mode 1, mode 2, and both');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
axis equal; grid on; view(3);
hold off;

%% =========================
% LOCAL FUNCTION (ported verbatim from Final_EI_STEP.m)
% =========================
function dirs = surfaceNormalsAtNodesFromTetMesh(mesh, nodeXYZ, nodeIDs)
% Estimate nodal outward surface normals from the free boundary of a
% tetrahedral mesh. Handles quadratic (10-node) tets by accumulating each
% boundary-face normal onto all of that face's nodes, including midside nodes,
% so midside surface nodes do not get a zero normal and a bad +Z fallback.
%
% Normals are area-weighted and flipped outward relative to the model centroid.
% For strongly non-convex geometry, visually verify the arrows.

E = mesh.Elements;

if size(E,1) < 4
    error("Expected tetrahedral elements with at least 4 nodes per element.");
end

nElem = size(E,2);
N = size(nodeXYZ,1);
centroid = mean(nodeXYZ,1);

% Corner-node face definitions (used to find boundary faces: a boundary face
% appears exactly once across all elements).
cornerFaceLocal = [
    1 2 3;
    1 2 4;
    2 3 4;
    1 3 4
];

% Full node sets per face for normal accumulation. Quadratic tets (>=10 local
% nodes) include midside nodes; assumed conventional local edge order:
%   5: edge 1-2, 6: edge 2-3, 7: edge 3-1,
%   8: edge 1-4, 9: edge 2-4, 10: edge 3-4.
if size(E,1) >= 10
    accumFaceLocal = {
        [1 2 3 5 6 7];      % face 1-2-3
        [1 2 4 5 9 8];      % face 1-2-4
        [2 3 4 6 10 9];     % face 2-3-4
        [1 3 4 7 10 8]      % face 1-3-4
    };
else
    accumFaceLocal = {
        [1 2 3];
        [1 2 4];
        [2 3 4];
        [1 3 4]
    };
end

% Build a table of all element faces keyed by sorted corner node IDs.
allFaceKeys = zeros(4*nElem, 3);
allFaceElem = zeros(4*nElem, 1);
allFaceNum  = zeros(4*nElem, 1);
row = 0;

for e = 1:nElem
    for f = 1:4
        row = row + 1;
        corners = E(cornerFaceLocal(f,:), e).';
        allFaceKeys(row,:) = sort(corners);
        allFaceElem(row) = e;
        allFaceNum(row) = f;
    end
end

[~, ~, ic] = unique(allFaceKeys, 'rows');
counts = accumarray(ic, 1);
isBoundaryRow = counts(ic) == 1;
boundaryRows = find(isBoundaryRow);

nodalNormals = zeros(N,3);
usedBoundaryNode = false(N,1);

for rr = boundaryRows(:).'
    e = allFaceElem(rr);
    f = allFaceNum(rr);

    cornerIDs = E(cornerFaceLocal(f,:), e).';
    p1 = nodeXYZ(cornerIDs(1),:);
    p2 = nodeXYZ(cornerIDs(2),:);
    p3 = nodeXYZ(cornerIDs(3),:);

    n = cross(p2 - p1, p3 - p1);   % area-weighted (length = 2*triangle area)
    if norm(n) < 1e-14
        continue;
    end

    faceNodeIDs = E(accumFaceLocal{f}, e).';
    faceNodeIDs = faceNodeIDs(faceNodeIDs >= 1 & faceNodeIDs <= N);
    faceNodeIDs = unique(faceNodeIDs(:));

    triCenter = mean(nodeXYZ(faceNodeIDs,:), 1);

    if dot(n, triCenter - centroid) < 0
        n = -n;   % flip outward
    end

    for a = 1:numel(faceNodeIDs)
        id = faceNodeIDs(a);
        nodalNormals(id,:) = nodalNormals(id,:) + n;
        usedBoundaryNode(id) = true;
    end
end

% Normalize the requested candidate-node normals.
dirs = nodalNormals(nodeIDs,:);
missing = false(numel(nodeIDs),1);

for i = 1:size(dirs,1)
    ni = norm(dirs(i,:));
    if ni < 1e-12
        missing(i) = true;
    else
        dirs(i,:) = dirs(i,:) / ni;
    end
end

% Robust fallback: copy the nearest valid boundary normal rather than a fixed
% global +Z, which is much safer on curved/angled surfaces.
if any(missing)
    goodBoundaryIDs = find(usedBoundaryNode & vecnorm(nodalNormals,2,2) > 1e-12);

    if isempty(goodBoundaryIDs)
        warning("No usable surface normals found anywhere. Falling back to +Z for %d nodes.", sum(missing));
        dirs(missing,:) = repmat([0 0 1], sum(missing), 1);
    else
        goodXYZ = nodeXYZ(goodBoundaryIDs,:);
        goodNormals = nodalNormals(goodBoundaryIDs,:);
        goodNormals = goodNormals ./ vecnorm(goodNormals,2,2);

        missingIdx = find(missing);
        for jj = 1:numel(missingIdx)
            i = missingIdx(jj);
            xyz = nodeXYZ(nodeIDs(i),:);
            [~, nearestLocal] = min(vecnorm(goodXYZ - xyz, 2, 2));
            dirs(i,:) = goodNormals(nearestLocal,:);
        end
        warning("Surface-normal fallback used nearest valid boundary normals for %d nodes.", sum(missing));
    end
end
end