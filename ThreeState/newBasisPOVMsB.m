function newPOVMs = newBasisPOVMsB(V)
%Takes in isometry from HV basis to \tilde{H}\tilde{V} basis
%Returns POVMs  for (tilde) H, V, R and L given V, n <= 1 subspace
%Assumes 2 + 1 dimensional Bob
arguments
    V (:, 2) double {isUnitaryNonSquare(V)}
end
numBasis = size(V, 1); 
newPOVMs = cell(numBasis, 1); 
for i = 1:numBasis
    v = [V(i, 1); V(i, 2)]; 
    gamma = conj(v)*transpose(v); 
    gammaherm = 1/2*(gamma+gamma');
    newPOVMs{i} = blkdiag(gammaherm, 0); 
end
end
function isUnitaryNonSquare(M)
if ~all(ismembertol(real(M'*M), eye(2)), 'all') || ~ismembertol(sum(imag(M'*M), 'all'), 0) 
    throw(MException('newBasisPOVMs:nonUnitary', ...
        'Matrix [ V | unused ] is not unitary.'));
end
end