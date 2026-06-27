function [keepIdx, hist] = energyOptimizationTechnique(Phi, massWeights, nKeep)
% energyOptimizationTechnique  Iterative sensor placement via EOT (Heo et al., 1997)
%
% The Energy Optimization Technique (EOT) is structurally identical to
% Kammer's Effective Independence Method (EIM), but replaces the Fisher
% Information Matrix with a Kinetic Energy matrix KE = Psi' * Psi, where
%
%   Psi = diag(sqrt(massWeights)) * Phi
%
% diag(sqrt(massWeights)) is the per-row Cholesky-like mass scaling that
% corresponds to the Cholesky factor U in the paper's derivation
% (M = L*U, Psi = U*Phi).
%
% At every iteration the EOT vector is computed as:
%
%   EOT_i = sum_j [ Psi_bar * Lambda^(-1/2) ]_ij ^2
%
% where Psi_bar is the current reduced projection and Lambda, eigenvectors
% are extracted from KE_bar = Psi_bar' * Psi_bar.
% The candidate with the smallest EOT score is removed, provided removal
% does not make the energy matrix rank-deficient.
%
% INPUTS
%   Phi         (nCandidates x nModes)  projected modal matrix (from
%               buildProjectedModalMatrix, already column-normalised)
%   massWeights (nCandidates x 1)       per-candidate mass proxy (scalar
%               per DOF).  Pass ones(size(Phi,1),1) to recover unweighted
%               behaviour.
%   nKeep       scalar  target number of sensors
%
% OUTPUTS
%   keepIdx     (nKeep x 1)  row indices into the original Phi/massWeights
%               arrays that form the optimal sensor set
%   hist        struct with fields:
%                 .numRemaining  number of active candidates before removal
%                 .minEOT        EOT score of the removed candidate

% --- Input validation --------------------------------------------------
if nargin < 3
    error('energyOptimizationTechnique requires Phi, massWeights, and nKeep.');
end

[nCand, nModes] = size(Phi);

if numel(massWeights) ~= nCand
    error('massWeights must have one entry per row of Phi (%d rows).', nCand);
end

if nKeep >= nCand
    warning('nKeep (%d) >= nCandidates (%d). Returning all candidates.', nKeep, nCand);
    keepIdx = (1:nCand)';
    hist.numRemaining = nCand;
    hist.minEOT = NaN;
    return;
end

if nKeep < nModes
    warning(['nKeep (%d) < nModes (%d). The resulting KE matrix will be rank-deficient. ' ...
             'Consider increasing nKeep or reducing nModesUse.'], nKeep, nModes);
end

% --- Initialise --------------------------------------------------------
active = (1:nCand)';   % indices of candidates still in the running

hist.numRemaining = [];
hist.minEOT       = [];

% Mass-scaling vector: sqrt of the diagonal of M restricted to candidates.
% This corresponds to the Cholesky factor U in Heo et al. eq. (4).
sqrtM = sqrt(massWeights(:));

% --- Iterative elimination ---------------------------------------------
while numel(active) > nKeep

    % Mass-scaled mode shape matrix for current active set
    PsiBar = sqrtM(active) .* Phi(active, :);   % (n_active x nModes)

    % Kinetic energy matrix KE_bar = PsiBar' * PsiBar   (nModes x nModes)
    KEbar = PsiBar' * PsiBar;

    % Eigendecomposition of KE_bar
    [V, Lam] = eig(KEbar);          % KEbar * V = V * Lam
    lambda = diag(Lam);             % eigenvalues

    % Guard against near-zero eigenvalues (rank deficiency).
    % Use a pseudo-inverse-style threshold.
    lambdaThresh = max(abs(lambda)) * 1e-10;
    safeInvSqrtLam = zeros(size(lambda));
    nonzero = abs(lambda) > lambdaThresh;
    safeInvSqrtLam(nonzero) = 1 ./ sqrt(abs(lambda(nonzero)));

    % Scaled eigenvector matrix: PsiBar * V * Lambda^(-1/2)
    % Size: (n_active x nModes)
    scaledEvec = PsiBar * V * diag(safeInvSqrtLam);

    % EOT score for each active candidate: row-wise sum of squares
    % Eq. (8) from Heo et al. 1997
    eotScores = sum(scaledEvec .^ 2, 2);   % (n_active x 1)

    % Record history
    hist.numRemaining(end+1, 1) = numel(active);
    hist.minEOT(end+1, 1)       = min(eotScores);

    % --- Rank-deficiency guard -----------------------------------------
    % Identify the candidate with the lowest EOT score
    [~, removeLocal] = min(eotScores);

    % Check that removing this candidate does not make KEbar rank-deficient
    testActive = active;
    testActive(removeLocal) = [];

    if numel(testActive) >= nModes
        testPsi = sqrtM(testActive) .* Phi(testActive, :);
        testKE  = testPsi' * testPsi;
        if rank(testKE) < nModes
            % Removing lowest-EOT candidate would drop rank; remove next-lowest
            [~, sortOrder] = sort(eotScores, 'ascend');
            removed = false;
            for k = 1:numel(sortOrder)
                candidate = sortOrder(k);
                trialActive = active;
                trialActive(candidate) = [];
                trialPsi = sqrtM(trialActive) .* Phi(trialActive, :);
                trialKE  = trialPsi' * trialPsi;
                if rank(trialKE) >= nModes
                    removeLocal = candidate;
                    removed = true;
                    break;
                end
            end
            if ~removed
                warning('Cannot remove any candidate without rank deficiency. Stopping early.');
                break;
            end
        end
    end

    active(removeLocal) = [];
end

keepIdx = active;
end
