function Phi = buildProjectedModalMatrix(R, nodeIDs, dirs, modeIDs)

if isprop(R.ModeShapes, 'ux')
    ux = R.ModeShapes.ux(nodeIDs, modeIDs);
    uy = R.ModeShapes.uy(nodeIDs, modeIDs);
    uz = R.ModeShapes.uz(nodeIDs, modeIDs);
elseif isprop(R.ModeShapes, 'x')
    ux = R.ModeShapes.x(nodeIDs, modeIDs);
    uy = R.ModeShapes.y(nodeIDs, modeIDs);
    uz = R.ModeShapes.z(nodeIDs, modeIDs);
else
    error('Unsupported ModeShapes format. Inspect R.ModeShapes.');
end

nNodes = numel(nodeIDs);
nModes = numel(modeIDs);

Phi = zeros(nNodes, nModes);

for i = 1:nNodes
    ni = dirs(i,:).';
    for m = 1:nModes
        u = [ux(i,m); uy(i,m); uz(i,m)];
        Phi(i,m) = ni.' * u;
    end
end

for m = 1:nModes
    nm = norm(Phi(:,m));
    if nm > 0
        Phi(:,m) = Phi(:,m) / nm;
    end
end
end