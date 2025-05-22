function mu = mu_beta_funcMat(N,F,t,epsilon)
    %N: Total number of rounds
    %F: Table of expectation values
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation
    
    Fflat = reshape(F,[],1);
    mu_init = ones(size(Fflat));
    tflat = reshape(t,[],1);
    
    for index=1:numel(Fflat)        
        %Calculate initial mu
        mu_init(index) = mu_beta_expec_floor(N,Fflat(index),tflat(index),epsilon);
    end

    %Reshape mu_init
    mu = reshape(mu_init,size(F));
end

function mu = mu_beta_expec_floor(N,F,t,epsilon)
    %N: Total number of rounds
    %F: expectation value for frequency
    %t: Acceptance parameter
    %epsilon: Epsilon used for parameter estimation

    %Lower range of the accepted parameters
    l1 = helper_minus(floor(N*(F-t)))-1;

    %Upper range of the accepted parameters
    l2 = floor(N*(F+t));
    
    %Function used for root finding
    fun = @(x) max(betainc(helper_01(1-(F-t-x)), N-l1, l1+1,"upper") , betainc(helper_01(1-(F+t+x)), helper_minus(N-l2), l2+1)) -epsilon;
    
    %Initial guesses
    x0 = 1/sqrt(N); % [0, 1-F-t];%muold(N,epsilon);
    
    %Mu given by root
    mu = fzero(fun,x0);
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