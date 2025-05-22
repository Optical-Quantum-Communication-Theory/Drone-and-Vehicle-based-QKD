function [Vd,cn,Vtot,W] = detectorDecompFlag(Mat,NB,perfect)
    %detectorDecompFlag calculates the decomposition of a detector matrix
    %describing a passive linear optics setup into a loss and detection
    %part
    %input:
    % Mat :             2Nout x Nin (lossy setup) or Nout x Nin (perfect setup) matrix 
    %                   describing the setup
    % NB (1) :          photon number cut-off
    % perfect (true):   parameter specifying if setup contains losses
    %                   (perfect = false) or not (perfect = true)
    %output
    % Vd :              Nout x Nin matrix describing the equivalent perfect 
    %                   detection setup
    % cn :              value needed for subspace estimation
    % Vtot :            blockdiagonal matrix containing Vd and Vl corresponding to the
    %                   isometries acting only on detection and loss
    % W :               isometry generating losses

    arguments
        Mat (:,:) double 
        NB (1,1) double {mustBeInteger,mustBePositive} = 1;
        perfect (1,1) logical = true;
    end
    
    % Perfect case with losses, i.e. Mat= Nout x Nin
    if perfect == true
        %number of input and output modes
        [Nout, Nin] =  size(Mat);
        
        %Iterate over all permutations of rows, corresponds to all possible
        %combinations of detectors
        vecPerms = flip(perms(1:Nout));
        
        %Store results in matrix
        SigmaMat = zeros(size(vecPerms,1),Nout/Nin);
        
        for index = 1:size(vecPerms,1)
            VMat = Mat(vecPerms(index,:),:);
            SigmaMat(index,:) = singularValFlag(VMat,Nin,Nout);
        end
    
        %Calculate subspace estimation for all possible combinations
        subspaceEstMat = sum(SigmaMat.^(2*(NB+1)),2);
        
        %Take minimum of subspace estimation matrix and make sure it is between
        %0 and 1
        cn = min(max(1 - min(subspaceEstMat),0),1);

        %For perfect case Vd = Mat, W = (id; 0), Vtot = Vd + 0
        Vd = Mat;
        W = [eye(Nin);zeros(Nin)];
        Vtot = blkdiag(Vd,zeros(Nout,Nin));

    % Imperfect case with losses, i.e. Mat= 2Nout x Nin
    else
        %number of input and output modes
        [Noutfull, Nin] =  size(Mat);
        
        %Outputmodes corresponding to a detection
        Nout = Noutfull/2;
        
        %Select parts of transformation resulting in detections only
        B = Mat(1:Nout,:);
        
        %Calculate norm of all vectors in B
        normB = vecnorm(B);
        
        %Normalise all vectors in B
        Vdinit = normalize_cols(B);
        
        %Perform svd for polar decomp
        [U1,S1,T1] = svd(Vdinit,"econ");
        
        %calculate polar decomp
        Vd = U1*T1';
        PB = T1*S1*T1';
        
        %Select parts of transformation resulting in losses only
        C = Mat(Nout+1:end,:);
        
        %Calculate norm of all vectors in B
        normC = vecnorm(C);
        
        %Normalise all vectors in B
        Vlinit = normalize_cols(C);
        
        %Perform svd for polar decomp
        [U2,S2,T2] = svd(Vlinit,'econ');
        
        %calculate polar decomp
        VC = U2*T2';
        PC = T2*S2*T2';
    
        %create total isometry acting on detections and losses seperately
        Vtot = blkdiag(Vd,VC);
    
        %Create isometry generating losses
        W = [PB*diag(normB); PC*diag(normC)];
    
        %Calculate bound for subspace estimation
        %number of rows in Vd
        numrows = size(Vd,1);
        
        %Iterate over all permutations of rows, corresponds to all possible
        %combinations of detectors
        vecPerms = flip(perms(1:numrows));
        
        %Store results in matrix
        SigmaMat = zeros(size(vecPerms,1),Nout/Nin);
        
        for index = 1:size(vecPerms,1)
            VMat = Vd(vecPerms(index,:),:);
            SigmaMat(index,:) = singularValFlag(VMat,Nin,Nout);
        end
        
        %Calculate subspace estimation for all possible combinations
        subspaceEstMat = sum(SigmaMat.^(2*(NB+1)),2);
        
        %Take minimum of subspace estimation matrix and make sure it is between
        %0 and 1
        cn = min(max(1 - min(subspaceEstMat),0),1);
    end
end

%Helper function calculating the subsapce estimation

function sigma = singularValFlag(VB,Nin,Nout)
    arguments
            VB (:,:) double 
            Nin (1,1) double {mustBeInteger,mustBePositive};
            Nout (1,1) double {mustBeInteger,mustBePositive};
    end
    %Calculate bound for flagstate squasher
    kmax = floor(Nout/Nin);
    
    %Store maximal singular values of each block in a vector
    %If Nout =/= kmax Nin, last block is not square 
    if kmax*Nin == Nout
        sigma = zeros(1,kmax);
    else
        sigma = zeros(1,kmax+1);
    end
    
    %Take SVD of all square blocks
    for indexk = 1:kmax
        Gblock = VB((indexk-1)*Nin+1:(indexk)*Nin,:);
        s = svd(Gblock);
        sigma(indexk) = max(s);
    end
    
    %Take SVD of last non-square block
    if kmax*Nin ~= Nout
        Gblocklast = VB(kmax*Nin+1:end,:);
        s = svd(Gblocklast);
        sigma(kmax+1) = max(s);
    end
end