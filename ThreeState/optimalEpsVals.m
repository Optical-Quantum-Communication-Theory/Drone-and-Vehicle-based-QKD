function epsilons = optimalEpsVals(N, n_sift, t_sift, dimA, epsSound, epsATfrac, opt)
%OPTIMALEPSVALS Quickly compute optimal epsilon values for finite size QKD
% Input:
% * N : Total Signal Count
% * n_sift : Generation signals that survive sifting
% * t_sift : Tolerance post sifting
% * dimA : Alice's system dimension
% * epsSound : Total epsilon security
% * epsATfrac : epsilonAT/epsilonSound. Fed as fraction for optimization purposes
% * ignoreExceed : Default 0. If 1, then any epsilon_i > epsilonSound will be allowed to be returned. 
% Output:
% * epsBar, epsEC, epsPA, epsAT
% Note epsSound = epsBar + epsEC + 1/2*(epsPA + epsAT)
% Default throws error if any of above epsilons exceeds total epsilon security.
% Aodhan Corrigan 31/10/2023
% Added output options 24/11/23
arguments
    N (1,1) double {mustBeInteger, mustBeNonnegative}
    n_sift (1,1) double {mustBeInteger, mustBeNonnegative}
    t_sift (1,1) double {mustBeNonnegative}
    dimA (1,1) double {mustBeInteger, mustBePositive}
    epsSound (1,1) double {mustBePositive}
    epsATfrac (1,1) double {mustBeInRange(epsATfrac, 0, 1)}
    opt.ignoreExceed (1,1) {mustBeInteger, mustBeInRange(opt.ignoreExceed, 0, 1)} = 0 %ignore if 1
    opt.verbosity (1,1) {mustBeInteger, mustBeInRange(opt.verbosity, 0, 2)} = 0
end

epsAT = epsSound*epsATfrac;

Delta = 2*sqrt(n_sift-N*t_sift)*log2(1+2*dimA);
r = 3 - 2*log2(epsAT);
phi = @(x) -2*log2(1-x) + Delta*sqrt(r - 2*log2(x)); %x=2*epsBar/epsAT
d = 1e-12; %Pertubation since phi(x) has vertical asymptotes at 0, 1.
x = fminbnd(phi, d, 1-d);

try
    epsBar = x*epsAT/2;
    epsEC = epsSound - epsAT;
    epsPA = epsAT - 2*epsBar;
    checkExceeds([epsBar, epsEC, epsPA], epsSound);
    epsilons =  [epsBar, epsEC, epsPA, epsAT];
catch error
    switch opt.ignoreExceed
        case 0
            rethrow(error)
        case 1
            warning('optimalEpsVals:ExceedSecurityWarn', ...
                'One or more security parameters exceeding total security. Ignoring error.')
            %Allow unphysical epsilons to continue with warning.
    end
end
switch opt.verbosity
    case 1
        fprintf('%s\n', 'Calculated optimal Epsilon values.')
    case 2
        fprintf('Optimal Epsilon weights epsBar epsEC epsPA epsAT: %.2f, %.2f, %.2f, %.2f \n', ...
            [epsBar/epsSound, epsEC/epsSound, epsPA/epsSound, epsAT/epsSound]); 
end
end

function checkExceeds(epsilon,epsSound)
%Verify none of the epsilon parameters exceeds the total key security.
%First checks constraint violated, then checks if due to numerical imprecision.
if any(epsilon > epsSound) && ~all(ismembertol(epsilon, epsSound))
    throw(MException('optimalEpsVals:ExceedSecurity', ...
        'One or more security parameters exceeding total security.'));
end
end

