function orbitGuessNew = haloUpdate(orbitHistory, continuationParameter)

n = numel(orbitHistory);

orbitGuessNew = copy(orbitHistory(n));

% Zero order approx
orbitGuessNew.x0(3) = continuationParameter;
% First order approx
if(n~=1)
    delta = orbitHistory(n).x0 - orbitHistory(n-1).x0;
    dxdp = delta([1 5]) / delta(3);
    dpNew = (orbitGuessNew.x0(3) - orbitHistory(n).x0(3));
    orbitGuessNew.x0([1 5]) = orbitGuessNew.x0([1 5]) + dxdp * dpNew;
    orbitGuessNew.x0([2 4 6]) = zeros(3,1); % force elements to 0
end