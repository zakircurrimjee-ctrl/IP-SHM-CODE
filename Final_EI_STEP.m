%% FAST EI-based Optimal Sensor Placement from CAD in MATLAB
% Requires: MATLAB + PDE Toolbox
% Surface-only EI version:
% - modal solve still uses the full volumetric mesh
% - EI is computed only on surface nodes from candidate faces
% - sensor measurement direction = local outward surface normal (see Section USER INPUTS)

clear; clc; close all;

%% =========================
% 0) SELECT CAD FILE
% ==========================

cadFile = "/Users/zakir/Desktop/Thesis/STEP/Test_EIGHT_Support_FUSED_MODULAR_3DIA_BetterENDS.STEP";

%% =========================
% UNIT HANDLING
% ==========================
cadLengthUnit = "m";

switch lower(cadLengthUnit)
    case "m"
        lengthScaleToMeters = 1.0;
    case "mm"
        lengthScaleToMeters = 1e-3;
    case "cm"
        lengthScaleToMeters = 1e-2;
    case "inch"
        lengthScaleToMeters = 0.0254;
    otherwise
        error('Unsupported cadLengthUnit. Use "m", "mm", "cm", or "inch".');
end

%% =========================
% USER INPUTS
% ==========================
fastMode = true;

    % Material properties - basic PLA
    E   = 2.66e9;    % Young's modulus, Pa
    nu  = 0.36;      % Poisson's ratio
    % EFFECTIVE (homogenized) density, kg/m^3.
    % The mesh treats the wall as solid continuum (nominal volume 605.44 cm^3),
    % but the as-printed wall is porous. rho here is the smeared value that
    % makes the SOLID model carry the part's REAL mass:
    %   rho_eff = as-printed mass / nominal solid volume
    %           = 750.91g / 605.44 cm^3 = 1240.27   (~60% of solid PLA 1250)
    % Voids are smeared UNIFORMLY: valid for global modes, does not capture
    % local mass variation if infill is spatially concentrated.
    % Mass 750.91 g is the measured (weighed) part mass.
    rho = 1240.27;       % effective density, kg/m^3 (was 1250 = solid PLA)

% % Material properties - tough PLA
% E   = 2.4e9;    % Young's modulus, Pa
% nu  = 0.35;      % Poisson's ratio
% rho = 1210;      % density, kg/m^3
% % --------------------------
% % FACE SELECTION OPTIONS
% % --------------------------
fixedFaces = [17 30];
candidateFaceMode = "manual";
candidateFacesManual = [1 26];;   % <-- your allowed sensor surface

forbiddenFaces = [];



% --------------------------
% SENSOR MEASUREMENT TYPE
% --------------------------
% Alternatives considered during development: local radial direction, the
% global X / Y / Z uniaxial directions, a user-defined custom axis, and an
% any-of-X/Y/Z uniaxial option. surfaceNormal was chosen because a uniaxial
% accelerometer glued to the curved tube wall measures along the local outward
% surface normal, so this is the physically correct measurement model here.
sensorMode = "surfaceNormal";

if fastMode
    nSensors = 2;
    freqRange = [-1, 2000]*2*pi;
    nModesUse = 2;

    minSpacing = 0.01;                 % 10 mm
    supportExclusionDistance = 0.01;    % 10 mm
    edgeExclusionDistance = 0.005;       % 5 mm

    amplitudeThresholdFrac = 0.05;
    candidateDownsampleStep = 1;
else
    nSensors = 2;
    freqRange = [-1, 2000]*2*pi;
    nModesUse = 2;

    minSpacing = 0.02;
    supportExclusionDistance = 0.01;
    edgeExclusionDistance = 0.005;

    amplitudeThresholdFrac = 0.10;
    candidateDownsampleStep = 1;
end

% --------------------------
% MESH ORDER + CONVERGENCE OPTIONS  (improvement #3)
% --------------------------
% Quadratic tetrahedral elements are recommended for bending-dominated modal
% analysis. Use "linear" only as a fallback if quadratic meshing is too costly.
meshGeometricOrder = "quadratic";

% useMeshConvergenceStudy = true:
%   Sweep HmaxCAD_list, track the first nModesConvergenceCheck natural
%   frequencies, and accept the first mesh whose frequencies change less than
%   meshConvergenceRelTol versus the next finer mesh. Writes a CSV table.
% useMeshConvergenceStudy = false:
%   Old fast behaviour: use the first mesh under the node target.
useMeshConvergenceStudy = true;

HminFraction          = 0.020;   % Hmin = HminFraction * Hmax
nodeCeiling           = 800000; % warn if a mesh exceeds this many nodes
meshConvergenceRelTol = 0.02;   % 2% relative frequency change = converged
nModesConvergenceCheck = max(3, nModesUse); % modes tracked in the study

%% =========================
% 1) IMPORT CAD + CREATE MODAL FEM MODEL
% ==========================
gm = fegeometry(cadFile);

% Scale CAD geometry to meters so mesh sizes, material properties, and
% frequencies are physically consistent.
if lengthScaleToMeters ~= 1
    gm = scale(gm, lengthScaleToMeters);
end

% --- Translate so the (minX, minY, minZ) bounding-box corner sits at origin ---
% Makes the whole part live in the positive octant, so reported sensor
% coordinates are positive and measurable from a single physical landmark.
bbMin = min(gm.Vertices, [], 1);   % [minX minY minZ] in meters, post-scale
originShift = -bbMin;              % translation that moves that corner to (0,0,0)
gm = translate(gm, originShift);
fprintf('Applied origin shift (m): [% .6f % .6f % .6f]\n', originShift);

model = femodel(AnalysisType="structuralModal", Geometry=gm);

figure;
pdegplot(model.Geometry, FaceLabels="on", EdgeLabels="on", FaceAlpha=0.4);
title("Inspect geometry face and edge labels");
axis equal; view(3);

allFaces = 1:model.Geometry.NumFaces;

if any(~ismember(fixedFaces, allFaces))
    error("One or more fixedFaces are invalid. Valid face IDs are 1:%d.", model.Geometry.NumFaces);
end

if any(~ismember(forbiddenFaces, allFaces))
    error("One or more forbiddenFaces are invalid. Valid face IDs are 1:%d.", model.Geometry.NumFaces);
end

switch lower(candidateFaceMode)
    case "all"
        candidateFaces = setdiff(allFaces, fixedFaces);
        candidateFaces = setdiff(candidateFaces, forbiddenFaces);

    case "manual"
         candidateFaces = unique(candidateFacesManual);

    forbiddenFaces = setdiff(allFaces, ...
                            [fixedFaces candidateFaces]);

    candidateFaces = setdiff(candidateFaces, fixedFaces);

    otherwise
        error('candidateFaceMode must be "all" or "manual".');
end

if isempty(candidateFaces)
    error("Candidate face list is empty after removing fixed/forbidden faces.");
end

fprintf("All faces: ");
disp(allFaces);

fprintf("Fixed faces: ");
disp(fixedFaces);

fprintf("Forbidden faces: ");
disp(forbiddenFaces);

fprintf("Candidate face mode: %s\n", candidateFaceMode);
fprintf("Candidate faces being used: ");
disp(candidateFaces);

model.MaterialProperties = materialProperties( ...
    YoungsModulus=E, ...
    PoissonsRatio=nu, ...
    MassDensity=rho);

for k = 1:numel(fixedFaces)
    model.FaceBC(fixedFaces(k)) = faceBC(Constraint="fixed");
end

%% =========================
% 2) MESH  (improvement #3: quadratic elements + convergence study)
% ==========================
% These values are in original CAD units.
% If cadLengthUnit = "mm", these are millimeters.
% Quadratic tets have many more DOFs than linear, so start coarser.
HmaxCAD_list = [10 5 3 2 1.2 0.9 0.6 0.3 0.1];

mesh = [];
lastGoodModel = [];
lastGoodMesh = [];
lastGoodH = NaN;
selectedHmaxCAD = NaN;

if useMeshConvergenceStudy
    fprintf('\n=== Mesh convergence study enabled ===\n');
    fprintf('Tracking the first %d modes. Relative tolerance = %.3g\n', ...
        nModesConvergenceCheck, meshConvergenceRelTol);
    fprintf('Mesh geometric order: %s\n', meshGeometricOrder);

    convFreqsHz = NaN(numel(HmaxCAD_list), nModesConvergenceCheck);
    convNodes   = NaN(numel(HmaxCAD_list), 1);
    convSuccess = false(numel(HmaxCAD_list), 1);

    for j = 1:numel(HmaxCAD_list)
        H = HmaxCAD_list(j) * lengthScaleToMeters;
        Hmin = H * HminFraction;

        fprintf('\nConvergence mesh %d/%d: HmaxCAD=%.4g %s, Hmax=%.6g m, order=%s\n', ...
            j, numel(HmaxCAD_list), HmaxCAD_list(j), cadLengthUnit, H, meshGeometricOrder);

        try
            testModel = generateMesh(model, Hmax=H, Hmin=Hmin, ...
                GeometricOrder=meshGeometricOrder);
            testMesh = testModel.Mesh;

            nNodes = size(testMesh.Nodes, 2);
            convNodes(j) = nNodes;
            convSuccess(j) = true;

            if nNodes > nodeCeiling
                warning(['Mesh has %d nodes (> %d). Modal solve may be slow ', ...
                    'or run out of memory. Increase Hmax if possible.'], nNodes, nodeCeiling);
            end

            % Lightweight modal solve just to read frequencies for the table.
            freqsHz_j = quickModalFrequenciesHz(testModel, freqRange, nModesConvergenceCheck);

            nCopy = min(numel(freqsHz_j), nModesConvergenceCheck);
            convFreqsHz(j,1:nCopy) = freqsHz_j(1:nCopy);

            fprintf('  Success: nodes=%d, first frequency=%.4f Hz\n', ...
                nNodes, convFreqsHz(j,1));

        catch ME
            fprintf('  Failed convergence solve: %s\n', ME.message);
            continue;
        end
    end

    validRows = find(convSuccess);
    if isempty(validRows)
        error(['All mesh attempts failed. Add finer entries to HmaxCAD_list ', ...
            'and double-check cadLengthUnit ("%s").'], cadLengthUnit);
    end

    % Build and print convergence table.
    T = table();
    T.HmaxCAD = HmaxCAD_list(:);
    T.Hmax_m  = HmaxCAD_list(:) * lengthScaleToMeters;
    T.Nodes   = convNodes;
    for m = 1:nModesConvergenceCheck
        T.(sprintf('f%d_Hz', m)) = convFreqsHz(:,m);
    end

    % Max relative frequency change (over nModesUse modes) vs the next finer mesh.
    maxRelChangeVsNext = NaN(numel(HmaxCAD_list), 1);
    for j = 1:numel(HmaxCAD_list)-1
        if convSuccess(j) && convSuccess(j+1)
            fCoarse = convFreqsHz(j,1:nModesUse);
            fFine   = convFreqsHz(j+1,1:nModesUse);
            ok = isfinite(fCoarse) & isfinite(fFine) & abs(fFine) > eps;
            if any(ok)
                maxRelChangeVsNext(j) = max(abs(fFine(ok) - fCoarse(ok)) ./ abs(fFine(ok)));
            end
        end
    end
    T.MaxRelChangeVsNext = maxRelChangeVsNext;

    disp('Mesh convergence table:');
    disp(T);
    writetable(T, 'mesh_convergence_frequencies.csv');

    convergedIdx = find(maxRelChangeVsNext <= meshConvergenceRelTol, 1, 'first');
    if isempty(convergedIdx)
        selIdx = validRows(end);   % finest successful mesh as fallback
        warning(['No mesh pair met the %.3g relative convergence tolerance. ', ...
            'Using the finest successful mesh: HmaxCAD=%.4g %s.'], ...
            meshConvergenceRelTol, HmaxCAD_list(selIdx), cadLengthUnit);
    else
        selIdx = convergedIdx + 1;  % finer mesh of the first converged pair
        fprintf(['Selected mesh from convergence: HmaxCAD=%.4g %s ', ...
            '(finer mesh of first converged pair).\n'], ...
            HmaxCAD_list(selIdx), cadLengthUnit);
    end

    selectedHmaxCAD = HmaxCAD_list(selIdx);
    H = selectedHmaxCAD * lengthScaleToMeters;
    Hmin = H * HminFraction;

    fprintf('\nGenerating final selected mesh: HmaxCAD=%.4g %s, Hmax=%.6g m, order=%s\n', ...
        selectedHmaxCAD, cadLengthUnit, H, meshGeometricOrder);

    lastGoodModel = generateMesh(model, Hmax=H, Hmin=Hmin, ...
        GeometricOrder=meshGeometricOrder);
    lastGoodMesh  = lastGoodModel.Mesh;
    lastGoodH     = H;

else
    % Fast path: first mesh that successfully generates at the chosen order.
    for j = 1:numel(HmaxCAD_list)
        H = HmaxCAD_list(j) * lengthScaleToMeters;
        Hmin = H * HminFraction;

        fprintf('\nTrying mesh: HmaxCAD=%.4g %s, Hmax=%.6g m, Hmin=%.6g m\n', ...
            HmaxCAD_list(j), cadLengthUnit, H, Hmin);

        try
            testModel = generateMesh(model, Hmax=H, Hmin=Hmin, ...
                GeometricOrder=meshGeometricOrder);
        catch ME
            fprintf('  Failed: %s\n', ME.message);
            continue;
        end

        testMesh = testModel.Mesh;
        nNodes = size(testMesh.Nodes, 2);
        fprintf('  Success: nodes=%d\n', nNodes);

        lastGoodModel = testModel;
        lastGoodMesh  = testMesh;
        lastGoodH     = H;
        selectedHmaxCAD = HmaxCAD_list(j);
        break;
    end
end

if isempty(lastGoodMesh)
    error("All mesh attempts failed. Your CAD may need repair/simplification, or try smaller HmaxCAD values.");
end

model = lastGoodModel;
mesh = lastGoodMesh;

allNodes = mesh.Nodes.';

modelMin = min(allNodes, [], 1);
modelMax = max(allNodes, [], 1);
modelSize = modelMax - modelMin;

fprintf('\nCAD unit selected: %s\n', cadLengthUnit);
fprintf('Final mesh nodes in full volumetric mesh: %d\n', size(allNodes,1));
fprintf('Geometry size after scaling: X=%.6g m, Y=%.6g m, Z=%.6g m\n', ...
    modelSize(1), modelSize(2), modelSize(3));
fprintf('Final successful Hmax: %.6g m (%s elements)\n', lastGoodH, meshGeometricOrder);
%% =========================
% 2b) VISUALISE MESH
% ==========================
% Three subplots in one figure:
%   Left   - full surface mesh wireframe (all nodes and edges)
%   Centre - nodes coloured by face ID (fixed / forbidden / candidate)
%   Right  - element quality histogram (aspect-ratio proxy via edge lengths)

meshNodes    = mesh.Nodes.';          % Nx3, all nodes
meshElements = mesh.Elements.';       % MxK, all volume elements (tetrahedral)

nMeshNodes    = size(meshNodes, 1);
nMeshElements = size(meshElements, 1);

fprintf('Mesh nodes: %d  |  Mesh elements: %d\n', nMeshNodes, nMeshElements);

% --- Collect surface nodes per face category ---
fixedNodes     = [];
forbiddenNodes = [];
candidateNodes = [];

for fID = allFaces
    faceNodeIDs = findNodes(mesh, 'region', 'Face', fID);
    if ismember(fID, fixedFaces)
        fixedNodes     = [fixedNodes,     faceNodeIDs]; %#ok<AGROW>
    elseif ismember(fID, forbiddenFaces)
        forbiddenNodes = [forbiddenNodes, faceNodeIDs]; %#ok<AGROW>
    else
        candidateNodes = [candidateNodes, faceNodeIDs]; %#ok<AGROW>
    end
end

fixedNodes     = unique(fixedNodes);
forbiddenNodes = unique(forbiddenNodes);
candidateNodes = unique(candidateNodes);

% --- Figure ---
figMesh = figure('Name', 'Mesh Visualisation', 'Color', 'w', ...
    'Position', [100 100 1000 480]);

% ---- Subplot 1: Surface wireframe ----
ax1 = subplot(1,2,1, 'Parent', figMesh);
pdeplot3D(model.Mesh, 'FaceAlpha', 0.15);
title(ax1, sprintf('Surface mesh\n%d nodes | %d elements', ...
    nMeshNodes, nMeshElements), 'FontSize', 10);
axis(ax1, 'equal'); view(ax1, 3); grid(ax1, 'on');
xlabel(ax1,'X (m)'); ylabel(ax1,'Y (m)'); zlabel(ax1,'Z (m)');

% ---- Subplot 2: Nodes coloured by face category ----
ax2 = subplot(1,2,2, 'Parent', figMesh);
hold(ax2, 'on');

if ~isempty(candidateNodes)
    c = meshNodes(candidateNodes, :);
    scatter3(ax2, c(:,1), c(:,2), c(:,3), 4, [0.2 0.7 0.3], 'filled', ...
        'DisplayName', 'Candidate faces');
end
if ~isempty(forbiddenNodes)
    f2 = meshNodes(forbiddenNodes, :);
    scatter3(ax2, f2(:,1), f2(:,2), f2(:,3), 10, [1 0.6 0], 'filled', ...
        'DisplayName', 'Forbidden faces');
end
if ~isempty(fixedNodes)
    f1 = meshNodes(fixedNodes, :);
    scatter3(ax2, f1(:,1), f1(:,2), f1(:,3), 18, [0.85 0.1 0.1], 'filled', ...
        'DisplayName', 'Fixed faces');
end

hold(ax2, 'off');
legend(ax2, 'Location', 'best', 'FontSize', 8);
title(ax2, 'Surface nodes by face category', 'FontSize', 10);
axis(ax2, 'equal'); view(ax2, 3); grid(ax2, 'on');
xlabel(ax2,'X (m)'); ylabel(ax2,'Y (m)'); zlabel(ax2,'Z (m)');



sgtitle(figMesh, 'Mesh inspection — verify before running EI', 'FontSize', 12, 'FontWeight', 'bold');

%% =========================
% 3) SOLVE MODAL PROBLEM
% ==========================
tic;
R = solve(model, FrequencyRange=freqRange);
fprintf('Modal solve time: %.2f s\n', toc);

freqs = R.NaturalFrequencies;
freqs_Hz = freqs / (2*pi);

if isempty(freqs)
    error("No modes found. Increase frequency range or check fixed face IDs.");
end

fprintf("Computed %d modes in the selected frequency range.\n", numel(freqs));

nShow = min(20, numel(freqs_Hz));
modeTable = table((1:nShow).', freqs_Hz(1:nShow), ...
    'VariableNames', {'ModeNumber','Frequency_Hz'});
disp(modeTable);

figure;
plot(1:numel(freqs_Hz), freqs_Hz, 'o-');
xlabel('Mode number');
ylabel('Frequency (Hz)');
title('Natural frequencies');
grid on;

if numel(freqs) < nModesUse
    warning("Only %d modes found. Reducing nModesUse to %d.", numel(freqs), numel(freqs));
    nModesUse = numel(freqs);
end

modeIDs = 1:nModesUse;

%% =========================
% 3b) VISUALIZE MODE SHAPES
% ==========================
plotModes = true;

if plotModes
    nModesToPlot = min(4, numel(freqs));
    deformationScale = 1.0;

    for m = 1:nModesToPlot
        ux = R.ModeShapes.ux(:,m);
        uy = R.ModeShapes.uy(:,m);
        uz = R.ModeShapes.uz(:,m);

        U = [ux, uy, uz];
        umag = sqrt(ux.^2 + uy.^2 + uz.^2);

        figure;
        pdeplot3D(mesh, ...
            "ColorMapData", umag, ...
            "Deformation", U, ...
            "DeformationScaleFactor", deformationScale);

        title(sprintf('Mode %d - %.2f Hz', m, freqs_Hz(m)));
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis equal;
        view(3);
        colorbar;
    end
end


%% =========================
% 4) EXTRACT SURFACE-ONLY CANDIDATE NODES
% ==========================
candidateNodeIDs = [];

for f = candidateFaces
    ids = findNodes(mesh, "region", "Face", f);
    candidateNodeIDs = [candidateNodeIDs; ids(:)]; %#ok<AGROW>
end

candidateNodeIDs = unique(candidateNodeIDs);

forbiddenNodeIDs = [];
for f = forbiddenFaces
    ids = findNodes(mesh, "region", "Face", f);
    forbiddenNodeIDs = [forbiddenNodeIDs; ids(:)]; %#ok<AGROW>
end
forbiddenNodeIDs = unique(forbiddenNodeIDs);

fixedNodeIDs = [];
for f = fixedFaces
    ids = findNodes(mesh, "region", "Face", f);
    fixedNodeIDs = [fixedNodeIDs; ids(:)]; %#ok<AGROW>
end
fixedNodeIDs = unique(fixedNodeIDs);

candidateNodeIDs = setdiff(candidateNodeIDs, forbiddenNodeIDs);
candidateNodeIDs = setdiff(candidateNodeIDs, fixedNodeIDs);

if isempty(candidateNodeIDs)
    error("No surface candidate nodes found. Check candidateFaces / forbiddenFaces / fixedFaces.");
end

fprintf("Initial SURFACE candidate nodes: %d\n", numel(candidateNodeIDs));

candidateNodeIDs = candidateNodeIDs(1:candidateDownsampleStep:end);
fprintf("Surface candidates after downsampling: %d\n", numel(candidateNodeIDs));

candidateXYZ = allNodes(candidateNodeIDs, :);

candidateNodeIDs_beforeFilter = candidateNodeIDs;
candidateXYZ_beforeFilter = candidateXYZ;

%% =========================
% 4b) SAVE WORKSPACE FOR DOWNSTREAM SCRIPTS
% ==========================
% Saved before EI filtering so the full surface candidate list is available
% to the MSE crack-placement and strike-mapping scripts.
fprintf("Saving workspace for MSE crack placement script...\n");
save('ei_workspace.mat', ...
    'R', ...                              % full modal result (mode shapes)
    'allNodes', ...                       % all mesh node coordinates
    'mesh', ...                           % mesh object (for visualisation)
    'model', ...                          % FEM model (for pdegplot)
    'freqs_Hz', ...                       % natural frequencies (Hz)
    'nModesUse', ...                      % number of modes used
    'modeIDs', ...                        % which mode indices were selected
    'candidateNodeIDs_beforeFilter', ...  % surface node IDs (pre-filter)
    'candidateXYZ_beforeFilter', ...      % surface node XYZ coords (pre-filter)
    'cadFile', ...                        % original CAD file path (for reference)
    'lengthScaleToMeters');               % unit scale factor (for reference)
fprintf("Saved: ei_workspace.mat\n\n");

%% =========================
% 5) BUILD MEASUREMENT MATRIX FOR SELECTED SENSOR TYPE
% ==========================

% Sensor measurement direction = local outward surface normal. This is the
% correct model for a uniaxial accelerometer glued to the curved tube wall.
% The sign of the normal does not affect EI: flipping a row of Phi leaves
% Phi'*Phi (and the leverage/EI scores) unchanged.
if lower(sensorMode) ~= "surfacenormal"
    error('This script supports sensorMode = "surfaceNormal" only.');
end

sensorDirs = surfaceNormalsAtNodesFromTetMesh(mesh, allNodes, candidateNodeIDs);
Phi = buildProjectedModalMatrix(R, candidateNodeIDs, sensorDirs, modeIDs);

sensorAxisLabel = repmat("surfaceNormal", numel(candidateNodeIDs), 1);

fprintf("Sensor mode: %s\n", sensorMode);
fprintf("Measurement candidates used by EI: %d\n", size(Phi,1));

if isempty(Phi) || size(Phi,1) ~= numel(candidateNodeIDs)
    error("Phi construction failed or does not match candidate count.");
end

%% =========================
% 6) FILTER BASED ON PHYSICS & PRACTICALITY
% ==========================

% 6a) Remove nodes too close to supports
supportNodeIDs = fixedNodeIDs;
supportXYZ = allNodes(supportNodeIDs, :);

if ~isempty(supportXYZ)
    keep = true(size(candidateXYZ,1),1);

    for i = 1:size(candidateXYZ,1)
        dmin = min(vecnorm(supportXYZ - candidateXYZ(i,:), 2, 2));
        keep(i) = dmin >= supportExclusionDistance;
    end

    candidateNodeIDs = candidateNodeIDs(keep);
    candidateXYZ = candidateXYZ(keep,:);
    sensorDirs = sensorDirs(keep,:);
    Phi = Phi(keep,:);
end

fprintf("Surface candidates after support filter: %d\n", size(Phi,1));

% 6b) Remove nodes too close to edges of candidate faces
edgeNodeIDs = [];

for f = candidateFaces
    faceNodeIDs = findNodes(mesh, "region", "Face", f);

    for otherFace = 1:model.Geometry.NumFaces
        if otherFace == f
            continue;
        end

        otherFaceNodeIDs = findNodes(mesh, "region", "Face", otherFace);
        sharedNodes = intersect(faceNodeIDs, otherFaceNodeIDs);

        edgeNodeIDs = [edgeNodeIDs; sharedNodes(:)]; %#ok<AGROW>
    end
end

edgeNodeIDs = unique(edgeNodeIDs);
edgeXYZ = allNodes(edgeNodeIDs, :);

if ~isempty(edgeXYZ)
    keep = true(size(candidateXYZ,1),1);

    for i = 1:size(candidateXYZ,1)
        dmin = min(vecnorm(edgeXYZ - candidateXYZ(i,:), 2, 2));
        keep(i) = dmin >= edgeExclusionDistance;
    end

    candidateNodeIDs = candidateNodeIDs(keep);
    candidateXYZ = candidateXYZ(keep,:);
    sensorDirs = sensorDirs(keep,:);
    Phi = Phi(keep,:);
end

if isempty(candidateNodeIDs)
    error("All candidates removed after edge filtering. Lower edgeExclusionDistance.");
end

fprintf("Surface candidates after edge filter: %d\n", size(Phi,1));

% 6c) Remove nodes with weak modal observability
modalRMS = sqrt(mean(Phi.^2, 2));
keep = modalRMS >= amplitudeThresholdFrac * max(modalRMS);

candidateNodeIDs = candidateNodeIDs(keep);
candidateXYZ = candidateXYZ(keep,:);
sensorDirs = sensorDirs(keep,:);
Phi = Phi(keep,:);

if isempty(candidateNodeIDs)
    error("All candidates removed after amplitude filtering. Lower amplitudeThresholdFrac.");
end

fprintf("Surface candidates after amplitude filter: %d\n", size(Phi,1));

% 6d) Enforce minimum spacing
modalRMS = sqrt(mean(Phi.^2, 2));
[~, order] = sort(modalRMS, "descend");

selected = false(size(candidateNodeIDs));
acceptedIDs = [];

for idx = order(:).'
    if isempty(acceptedIDs)
        acceptedIDs(end+1) = idx; %#ok<SAGROW>
        selected(idx) = true;
    else
        d = vecnorm(candidateXYZ(acceptedIDs,:) - candidateXYZ(idx,:), 2, 2);
        if all(d >= minSpacing)
            acceptedIDs(end+1) = idx; %#ok<SAGROW>
            selected(idx) = true;
        end
    end
end

candidateNodeIDs = candidateNodeIDs(selected);
candidateXYZ = candidateXYZ(selected,:);
sensorDirs = sensorDirs(selected,:);
Phi = Phi(selected,:);

fprintf("Surface candidates after spacing filter: %d\n", size(Phi,1));

if size(Phi,1) < nSensors
    warning("Filtered candidate pool (%d) is smaller than requested number of sensors (%d). Reducing nSensors.", ...
        size(Phi,1), nSensors);
    nSensors = size(Phi,1);
end

fprintf("Candidates before filtering plot set: %d\n", size(candidateXYZ_beforeFilter,1));
fprintf("Filtered candidate pool used by EI: %d\n", size(candidateXYZ,1));


%% =========================
% 7) APPLY EI ON SURFACE NODES ONLY
% ==========================
tic;

% Best sensors: final EI survivors
[bestIdxRaw, eiHistory] = effectiveIndependence(Phi, nSensors);

bestIdxLocal = bestIdxRaw;

% Worst / medium sensors: selected by EI score over the feasible candidate pool
Aall = Phi;
Fall = Aall' * Aall;

if rcond(Fall) < 1e-12
    FallInv = pinv(Fall);
else
    FallInv = inv(Fall);
end

EIall = diag(Aall * FallInv * Aall');

% Sort candidates by EI score
% Low EI = poor / redundant candidate
% High EI = informative candidate
[~, orderAsc] = sort(EIall, "ascend");

worstIdxRaw = orderAsc;

midCenter = round(numel(orderAsc)/2);
midWindow = max(10*nSensors, nSensors);
midStart = max(1, midCenter - floor(midWindow/2));
midEnd = min(numel(orderAsc), midStart + midWindow - 1);
mediumIdxRaw = orderAsc(midStart:midEnd);

mediumIdxLocal = mediumIdxRaw(1:min(nSensors,numel(mediumIdxRaw)));
worstIdxLocal = worstIdxRaw(1:min(nSensors,numel(worstIdxRaw)));

fprintf('EI time: %.2f s\n', toc);

bestNodeIDs = candidateNodeIDs(bestIdxLocal);
bestXYZ = candidateXYZ(bestIdxLocal,:);
bestSensorDirs = sensorDirs(bestIdxLocal,:);
bestEIscore = EIall(bestIdxLocal);

mediumNodeIDs = candidateNodeIDs(mediumIdxLocal);
mediumXYZ = candidateXYZ(mediumIdxLocal,:);
mediumSensorDirs = sensorDirs(mediumIdxLocal,:);
mediumEIscore = EIall(mediumIdxLocal);

worstNodeIDs = candidateNodeIDs(worstIdxLocal);
worstXYZ = candidateXYZ(worstIdxLocal,:);
worstSensorDirs = sensorDirs(worstIdxLocal,:);
worstEIscore = EIall(worstIdxLocal);

fprintf("Best EI sensors: %d\n", size(bestXYZ,1));
fprintf("Medium EI sensors: %d\n", size(mediumXYZ,1));
fprintf("Worst EI sensors: %d\n", size(worstXYZ,1));

fprintf("\nMean EI scores:\n");
fprintf("Best:   %.6g\n", mean(bestEIscore));
fprintf("Medium: %.6g\n", mean(mediumEIscore));
fprintf("Worst:  %.6g\n", mean(worstEIscore));

%% =========================
% 7b) SENSOR-SET DIAGNOSTICS  (improvement #4)
% ==========================
% Quantify how well each chosen subset resolves the target modes:
%   det(Phi'Phi)  -> D-optimality (higher = better)
%   cond(Phi'Phi) -> conditioning (lower = better)
%   singular values, per-mode signal norms
%   auto-MAC      -> mode confusability (low off-diagonal = better)
% Note: Phi columns were normalized per mode in buildProjectedModalMatrix, so
% these are relative comparisons between the three sets, not absolute energies.
PhiBest   = Phi(bestIdxLocal,   :);
PhiMedium = Phi(mediumIdxLocal, :);
PhiWorst  = Phi(worstIdxLocal,  :);

infoBest   = sensorSetDiagnostics(PhiBest,   "Best EI");
infoMedium = sensorSetDiagnostics(PhiMedium, "Medium EI");
infoWorst  = sensorSetDiagnostics(PhiWorst,  "Worst EI");

%% =========================
% 8) MAP TO REAL SENSOR POSITIONS
% ==========================

bestAxisLabel = directionLabelsFromDirs(bestSensorDirs);
mediumAxisLabel = directionLabelsFromDirs(mediumSensorDirs);
worstAxisLabel = directionLabelsFromDirs(worstSensorDirs);

bestSensorTable = table( ...
    (1:numel(bestNodeIDs)).', ...
    bestNodeIDs(:), ...
    bestAxisLabel(:), ...
    bestXYZ(:,1), bestXYZ(:,2), bestXYZ(:,3), ...
    bestSensorDirs(:,1), bestSensorDirs(:,2), bestSensorDirs(:,3), ...
    bestEIscore(:), ...
    'VariableNames', {'SensorID','NodeID','Axis','X','Y','Z','nx','ny','nz','EIscore'});


worstSensorTable = table( ...
    (1:numel(worstNodeIDs)).', ...
    worstNodeIDs(:), ...
    worstAxisLabel(:), ...
    worstXYZ(:,1), worstXYZ(:,2), worstXYZ(:,3), ...
    worstSensorDirs(:,1), worstSensorDirs(:,2), worstSensorDirs(:,3), ...
    worstEIscore(:), ...
    'VariableNames', {'SensorID','NodeID','Axis','X','Y','Z','nx','ny','nz','EIscore'});

disp("Best EI sensor positions:");
disp(bestSensorTable);


disp("Worst EI sensor positions:");
disp(worstSensorTable);

writetable(bestSensorTable, "best_sensor_positions.csv");
writetable(worstSensorTable, "worst_sensor_positions.csv");

%% =========================
% 8b) SAVE EI RESULTS FOR MAC VALIDATION
% ==========================
% Saves the filtered mode shape matrix and best sensor indices so that
% compute_mac_validation.m can load them without re-running the FEM solve.

Phi_EI       = Phi;           % filtered projected mode shape matrix (nCandidates x nModes)
bestIdx_EI   = bestIdxLocal;  % indices into Phi_EI for the best sensors
bestNodeIDs_EI = bestNodeIDs; % physical node IDs of best sensors
bestXYZ_EI   = bestXYZ;       % XYZ coordinates of best sensors
freqs_Hz_EI  = freqs_Hz;      % natural frequencies (Hz)
nModesUse_EI = nModesUse;     % number of target modes

save('ei_osp_results.mat', ...
    'Phi_EI', ...
    'bestIdx_EI', ...
    'bestNodeIDs_EI', ...
    'bestXYZ_EI', ...
    'freqs_Hz_EI', ...
    'nModesUse_EI');

fprintf('\nSaved EI results to ei_osp_results.mat\n');
fprintf('  Phi_EI size:     %d candidates x %d modes\n', size(Phi_EI,1), size(Phi_EI,2));
fprintf('  Best sensors:    %d\n', numel(bestIdx_EI));

%% =========================
% 9) VISUALIZE SURFACE CANDIDATES + BEST / MEDIUM / WORST SENSORS
% ==========================
figure;
pdegplot(model.Geometry, FaceAlpha=0.10);
hold on;

hBefore = scatter3(candidateXYZ_beforeFilter(:,1), candidateXYZ_beforeFilter(:,2), candidateXYZ_beforeFilter(:,3), ...
    18, [0.6 0.6 1.0], 'filled', ...
    'MarkerFaceAlpha', 0.20, 'MarkerEdgeAlpha', 0.20);

hFiltered = scatter3(candidateXYZ(:,1), candidateXYZ(:,2), candidateXYZ(:,3), ...
    30, [0.0 0.8 0.0], 'filled', ...
    'MarkerFaceAlpha', 0.60, 'MarkerEdgeAlpha', 0.60);

hBest = scatter3(bestXYZ(:,1), bestXYZ(:,2), bestXYZ(:,3), ...
    120, 'r', 'filled');


hWorst = scatter3(worstXYZ(:,1), worstXYZ(:,2), worstXYZ(:,3), ...
     120, 'k', 'filled');

quiver3(bestXYZ(:,1), bestXYZ(:,2), bestXYZ(:,3), ...
        bestSensorDirs(:,1), bestSensorDirs(:,2), bestSensorDirs(:,3), ...
        0.03, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

 quiver3(mediumXYZ(:,1), mediumXYZ(:,2), mediumXYZ(:,3), ...
         mediumSensorDirs(:,1), mediumSensorDirs(:,2), mediumSensorDirs(:,3), ...
         0.03, 'm', 'LineWidth', 1.5, 'HandleVisibility', 'off');

 quiver3(worstXYZ(:,1), worstXYZ(:,2), worstXYZ(:,3), ...
         worstSensorDirs(:,1), worstSensorDirs(:,2), worstSensorDirs(:,3), ...
        0.03, 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --------------------------
% Add sensor ID labels
% --------------------------
labelOffset = 0.003;   % meters; increase if labels overlap

for i = 1:size(bestXYZ,1)
    text(bestXYZ(i,1)+labelOffset, bestXYZ(i,2)+labelOffset, bestXYZ(i,3)+labelOffset, ...
        sprintf('B%d', i), ...
        'Color', 'r', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
end
 
 for i = 1:size(worstXYZ,1)
     text(worstXYZ(i,1)+labelOffset, worstXYZ(i,2)+labelOffset, worstXYZ(i,3)+labelOffset, ...
         sprintf('W%d', i), ...
         'Color', 'k', ...
         'FontSize', 10, ...
         'FontWeight', 'bold', ...
         'BackgroundColor', 'w', ...
         'Margin', 1);
 end

legend([hBefore, hFiltered, hBest, hWorst], ...
       {"Surface candidates before filter", ...
        "Filtered candidates", ...
        "Best EI sensors", ...
        "Worst EI sensors"}, ...
       "Location", "best");

title(sprintf("Best  / Worst EI sensor positions, sensorMode = %s", sensorMode));
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);
%% =========================
% 10) SHOW EI ELIMINATION HISTORY
% ==========================
figure;
plot(eiHistory.numRemaining, eiHistory.minEI, '-o');
xlabel("Number of remaining candidates");
ylabel("Minimum EI value removed at each step");
title("EI elimination history");
grid on;



function labels = directionLabelsFromDirs(dirs)

labels = strings(size(dirs,1),1);

for i = 1:size(dirs,1)
    d = dirs(i,:);

    if norm(d - [1 0 0]) < 1e-8
        labels(i) = "x";
    elseif norm(d - [0 1 0]) < 1e-8
        labels(i) = "y";
    elseif norm(d - [0 0 1]) < 1e-8
        labels(i) = "z";
    elseif norm(d + [1 0 0]) < 1e-8
        labels(i) = "-x";
    elseif norm(d + [0 1 0]) < 1e-8
        labels(i) = "-y";
    elseif norm(d + [0 0 1]) < 1e-8
        labels(i) = "-z";
    else
        labels(i) = "custom/radial";
    end
end

end

%% ============================================================
% ADDED HELPER FUNCTIONS (improvements #1, #3, #4)
% ============================================================

function freqsHz = quickModalFrequenciesHz(model, freqRange, nModesWanted)
% Lightweight modal solve used ONLY by the mesh-convergence study.
% Whole-face fixed BCs are already applied to `model`.
%
% Speed: the study only needs the first few frequencies to judge convergence,
% so we ask the solver for the smallest nModesWanted modes instead of every
% mode up to freqRange(2). Computing N modes is far cheaper than a full-range
% solve on a fine mesh. Falls back to the range solve on older toolboxes.
try
    R = solve(model, ModalResults=nModesWanted);  %#ok<*NASGU>  % newer syntax
catch
    try
        R = solve(model, NumModes=nModesWanted);            % alternate newer syntax
    catch
        R = solve(model, FrequencyRange=freqRange);         % robust fallback
    end
end

omega = R.NaturalFrequencies;
omega = omega(isfinite(omega) & omega >= freqRange(1));
omega = sort(omega, "ascend");

freqsHz = omega / (2*pi);

if numel(freqsHz) > nModesWanted
    freqsHz = freqsHz(1:nModesWanted);
end
end


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


function info = sensorSetDiagnostics(PhiSub, setName)
% Quantify how well a chosen sensor subset resolves the target modes.
% PhiSub is (nSensors x nModes): rows = selected sensors, cols = target modes.

nModes = size(PhiSub, 2);

F = PhiSub.' * PhiSub;          % Fisher information matrix (nModes x nModes)

detF  = det(F);                 % D-optimality (EI maximizes this)
condF = cond(F);                % conditioning

sv    = svd(PhiSub);
minSV = min(sv);
maxSV = max(sv);

colNorms = vecnorm(PhiSub, 2, 1);  % per-mode signal magnitude at these sensors

% Auto-MAC of the mode-shape partition at these sensors.
MAC = zeros(nModes);
for i = 1:nModes
    for j = 1:nModes
        num = (PhiSub(:,i).' * PhiSub(:,j))^2;
        den = (PhiSub(:,i).' * PhiSub(:,i)) * (PhiSub(:,j).' * PhiSub(:,j));
        if den > 0
            MAC(i,j) = num / den;
        else
            MAC(i,j) = NaN;
        end
    end
end

offDiag   = MAC - diag(diag(MAC));
maxOffMAC = max(offDiag(:));

info.name      = setName;
info.F         = F;
info.detF      = detF;
info.condF     = condF;
info.sv        = sv;
info.minSV     = minSV;
info.maxSV     = maxSV;
info.colNorms  = colNorms;
info.MAC       = MAC;
info.maxOffMAC = maxOffMAC;

fprintf('\n--- Sensor-set diagnostics: %s ---\n', setName);
fprintf('  det(Phi''*Phi)  (D-opt, higher = better)   : %.6g\n', detF);
fprintf('  cond(Phi''*Phi) (lower = better)            : %.6g\n', condF);
fprintf('  min / max singular value of Phi_sub        : %.4g / %.4g\n', minSV, maxSV);
fprintf('  per-mode signal norm at sensors            : %s\n', mat2str(colNorms, 3));
fprintf('  max off-diagonal auto-MAC (lower = better) : %.4f\n', maxOffMAC);
fprintf('  auto-MAC matrix:\n');
disp(MAC);
end