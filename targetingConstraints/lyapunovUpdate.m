function orbitGuessNew = lyapunovUpdate(orbitHistory, continuationParameter)
n = numel(orbitHistory);

orbitGuessNew = copy(orbitHistory(n));

% Zero order approx
orbitGuessNew.x0(1) = continuationParameter;
% First order approx
if(n~=1)
    delta = orbitHistory(n).x0 - orbitHistory(n-1).x0;
    dxdp = delta(2:end) / delta(1);
    dpNew = (orbitGuessNew.x0(1) - orbitHistory(n).x0(1));
    orbitGuessNew.x0(2:end) = orbitGuessNew.x0(2:end) + dxdp * dpNew;
    orbitGuessNew.x0([2 3 4 6]) = zeros(4,1); % force elements to 0
end