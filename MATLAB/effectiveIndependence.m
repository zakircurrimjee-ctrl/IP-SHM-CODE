function [bestIdx, hist, mediumIdx, worstIdx] = effectiveIndependence(Phi, nKeep)

active = (1:size(Phi,1))';

hist.numRemaining = [];
hist.minEI = [];
hist.removedIdx = [];
hist.removedEI = [];

while numel(active) > nKeep
    A = Phi(active,:);
    F = A' * A;

    if rcond(F) < 1e-12
        Finv = pinv(F);
    else
        Finv = inv(F);
    end

    E = diag(A * Finv * A');

    [minVal, idxLocal] = min(E);

    removedGlobalIdx = active(idxLocal);

    hist.numRemaining(end+1,1) = numel(active);
    hist.minEI(end+1,1) = minVal;
    hist.removedIdx(end+1,1) = removedGlobalIdx;
    hist.removedEI(end+1,1) = minVal;

    active(idxLocal) = [];
end

% Best sensors are the final remaining candidates
bestIdx = active;

% Worst sensors are the first removed candidates
nWorst = min(nKeep, numel(hist.removedIdx));
worstIdx = hist.removedIdx(1:nWorst);

% Medium sensors are taken from the middle of the removal history
nRemoved = numel(hist.removedIdx);
nMedium = min(nKeep, nRemoved);

midStart = max(1, round(nRemoved/2 - nMedium/2));
midEnd = min(nRemoved, midStart + nMedium - 1);

mediumIdx = hist.removedIdx(midStart:midEnd);

end