function [muLower,muUpper] = mubetaUpperLower(N,F,t,logEpsilon)
    arguments        
        N (1,1) uint64 {mustBePositive}
        F (:,:,:) double {mustBeNonnegative, mustBeLessThanOrEqual(F,1)}
        t (:,:,:) double {mustBeNonnegative, mustBeLessThanOrEqual(t,1)}
        logEpsilon (1,1) double {mustBeNegative}
    end
    %N: Total number of rounds
    %F: Table of expectation values
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation
    
    %Store shape of input array
    shape = size(F);
    
    %Flatten arrays
    Fflat = F(:);
    muLower = ones(size(Fflat));
    muUpper = ones(size(Fflat));
    tflat = t(:);
    
    if logEpsilon >= -13
        %Covert log of epsilon back into epsilon
        epsilon = 10^(logEpsilon);

        for index=1:numel(Fflat)        
            %Calculate initial mu
            [muLower(index), muUpper(index)] = mu_beta_expec_floor(N,Fflat(index),tflat(index),epsilon);
        end
    
        %Reshape mus
        muLower = reshape(muLower,shape);
        muUpper = reshape(muUpper,shape);

    else
        %Implement bounds using a Gaussian bound on binomial distribution
        %from A Complete Proof of Universal Inequalities for the Distribution
        % Function of the Binomial Law 
        % https://doi.org/10.1137/S0040585X97986138

        for index=1:numel(Fflat)        
            %Calculate initial mu
            [muLower(index), muUpper(index)] = mu_Gauss(N,Fflat(index),tflat(index),logEpsilon);
        end
    
        %Reshape mus
        muLower = reshape(muLower,shape);
        muUpper = reshape(muUpper,shape);
    end
end

function [muLower,muUpper] = mu_Gauss(N,F,t,logEpsilon)
    %N: Total number of rounds  
    %F: expectation value for frequency
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation

    %Store symbolic version of epsilon
    epsSym = 10^(sym(logEpsilon));

    %Calculate inverse of normal distribution of epsilon
    phiInvEps = -sqrt(2)*erfcinv(2*epsSym);
    phiInv1MinusEps = -sqrt(2)*erfcinv(2*(1-epsSym));

    %Convert N into double
    Ndbl = double(N);

    %% Upper bound

    %Upper range of the accepted parameters
    kUpper = double(floor(N*(F+t)));
    cUpper = double(phiInvEps^2/2);

    funUpper = @(x) ((kUpper+1)*(log((kUpper+1)/Ndbl) - log(x)) + ...
                (Ndbl - kUpper - 1/Ndbl)*(log(1-(kUpper+1)/Ndbl) - log(1-x)) - cUpper);
    

    %% Lower bound

    %Lower range of the accepted parameters
    kLower = double(ceil(N*(F-t))-1);
    cLower = double(phiInv1MinusEps^2/2);
    
    %Functions used for root finding
    funLower= @(x) (kLower*(log(kLower/Ndbl) - log(x)) + ...
                (Ndbl - kLower)*(log(1-kLower/Ndbl) - log(1-x)) - cLower);
    
    
    %% Find bound by root finding
    cutoff = 1e-16;
    %Mu given by root
    if kLower <= 0
        muLower = min(max(double(F - t - 1 + nthroot(1-epsSym,sym(N))),0), abs(F-t));
    elseif kLower >= N
        muLower = double(F - t - nthroot(epsSym,sym(N)));
    else
        muLower = abs((F - t) - fzero(funLower,[eps(0),kLower/Ndbl]));
    end

    if kUpper >= N
        muUpper = min(max(double(nthroot(1-epsSym,sym(N)) - (F + t) ),0), abs(1 - (F+t))) ;
    elseif kUpper <= 0
        muUpper = double( 1 - nthroot(epsSym,sym(N)) - (F + t) );
    else
        muUpper = abs(fzero(funUpper,[kUpper/Ndbl,1-kUpper/Ndbl]) - (F + t));
    end
end


function [muLower,muUpper] = mu_beta_expec_floor(N,F,t,epsilon)
    %N: Total number of rounds  
    %F: expectation value for frequency
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation

    %Lower range of the accepted parameters
    l1 = helper_minus(ceil(N*(F-t)))-1;

    %Upper range of the accepted parameters
    l2 = double(floor(N*(F+t)));

    %Convert N into double
    Ndbl = double(N);
    
    %Function used for root finding
    funLower = @(x) betainc(helper_01(1-(F-t-x)), Ndbl-l1, l1+1,"upper") - epsilon;
    funUpper = @(x) betainc(helper_01(1-(F+t+x)), helper_minus(Ndbl-l2), l2+1) - epsilon;
    
    %Initial guess
    x0 = 1/sqrt(Ndbl);
    
    %Mu given by root
    muLower = max(fzero(funLower,x0),0);
    muUpper = fzero(funUpper,x0);
end

function value = helper_minus(x)
    %Function converts array to an array with non-negative entries
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        else
            value(index) = x(index);
        end
    end
end

function value = helper_01(x)
    %Function converts array to an array with entries between 0 and 1
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        elseif x(index) >1
            value(index) = 1;
        else
            value(index) = x(index);
        end
    end
end