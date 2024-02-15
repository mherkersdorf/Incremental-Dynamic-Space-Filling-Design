function [u, uProxy, yProxy, levels, visitsLevels] = IDSFIDARX(N, dim, T, Ts, lambda,  M, amplitudeConstraints, velocityConstraints, additionalConstraints, existingData, Lmax, n, gains)
% It follows an implementation of the Incremental Dynamic Space-Filling Design
% (IDS-FID) approach, a Design of Experiments technique for nonlinear dynamic
% systems. It aimes at achieving a space-filling distribution within the regressor
% space of a (proxy) model while simultaneously offering control over the excited
% spectrum by adjusting the hyperparameter lambda. Here, the proxy is established
% using a first-order ARX model.
%
% Author / Copyright: Oliver Nelles and Max Herkersdorf
%
% E-Mail: Max.Herkersdorf@uni-siegen.de     Oliver.Nelles@uni-siegen.de
%
% Date: 2024/22/01
arguments
    N (1,1) double                                                           % Number of data points.
    dim (1,1) double                                                         % Dimensions.
    T (1,:) double                                                           % Time constants from each input to the output.
    Ts (1,1) double                                                          % Sampling time in seconds.
    lambda (1,:) double                                                      % Value of the hyperparameter lambda for each input.
    M (1,:) double = ceil(N./(2*dim*T/Ts))                                   % Number of visitable amplitude levels.
    amplitudeConstraints (2,:) double = [0; 1]*ones(1,dim)                   % Minimum (1st row) and maximum (2nd row) amplitude value vor each input.
    velocityConstraints (2,:) double = [-(amplitudeConstraints(2,:) - ...
        amplitudeConstraints(1,:)); (amplitudeConstraints(2,:) - ...
        amplitudeConstraints(1,:))]                                          % Constraints on the rate of change for each input during decrease (1st row) and increase (2nd row). By default no constraints.
    additionalConstraints (1,:) string = 'no'                                % Further user-defined constraints; e.g., prohibit certain regions in the (proxy) regressr space, ban on combinations of amplitude levels ...
    existingData (:,:) double = amplitudeConstraints(1,:)                    % Implementation of existing data. If no existing data is on hand, point (0, 0, 0 ...) (in proxy regressor space) is employed to start the optimization.
    Lmax (1,:) double = ceil(dim*4.*T/Ts)                                    % Maximum sequence length for each input.
    n (1,:) double = 4*ones(1,dim)                                           % Exponent of the penalty term for each input.
    gains (1,:) double = ones(1,dim) / dim                                   % Gains from each input to the output. Note: Due to normalization gains need not necessarily reflect the real process gains.  Instead, they indicate the dominant regressors during optimization. Default gains equal 1/dim, ensuring that all inputs carry equal importance.
end

%% Initialization
% Prepare data for optimization.
assert(size(T,2) == dim && size(lambda,2) == dim && size(M,2) == dim && size(Lmax,2) == dim && size(gains,2) == dim && size(velocityConstraints,2) == dim && size(existingData,2) == dim && size(amplitudeConstraints,2) == dim, "Number of values per parameter must match the number of dimensions.")
[levels, J2, visitsLevels] = deal(cell(1,dim));
for jj = 1:dim
    % Create visitable amplitude levels.
    levels{1, jj} = linspace(0, 1, M(jj))';
    % Create penalty term.
    J2{1,jj} = (1:Lmax(jj)).^n(jj)';
    J2{1, jj} = normalize(J2{1, jj}, 'range');
    % Create visit counter for amplitude levels to incorporate
    % Latin hypercube-like approach.
     visitsLevels{1,jj} = deal(zeros(M(jj), 1));
    % Adjust velocity constraints due to normalization in range [0,1]
    velocityConstraints(1,jj)  = (-1)*rescale((-1)*velocityConstraints(1,jj),"InputMin",0,"InputMax",amplitudeConstraints(2,jj)-amplitudeConstraints(1,jj));
    velocityConstraints(2,jj)  = rescale(velocityConstraints(2,jj),"InputMin",0,"InputMax",amplitudeConstraints(2,jj)-amplitudeConstraints(1,jj));
end
% Normalize existing data between [0, 1] with regard to the possible maximum / minimum values in this dimension.
existingData = rescale(existingData,"InputMin",amplitudeConstraints(1,:),"InputMax",amplitudeConstraints(2,:));
uOpt = mat2cell(existingData, size(existingData, 1), repelem(1, dim));
% Create proxy model.
polesSys = exp(-Ts./T);                 gainsSys = (1 - polesSys).*gains;
polesSys = num2cell(polesSys);          zerosSys = num2cell(zeros(1,dim));
Sys = zpk(zerosSys, polesSys, gainsSys, Ts);

%% Optimization
iterations = ones(1,dim);
while any(cell2mat(cellfun(@length, uOpt, 'UniformOutput', false)) < N)                       % Optimize till preset signal length is reached for all inputs.
    [~, idx] = min(cell2mat(cellfun(@length, uOpt, 'UniformOutput', false)));                 % Chose input with the fewest  data points for optimization.
    % Number of already optimized data points per dimension.
    numberOptDataPoints = cell2mat(cellfun(@length, uOpt, 'UniformOutput', false));

    % Available sequences in this iteration for the current input.
    availableSequencesTemp = repmat(levels{idx}', Lmax(idx),1);

    % Exclude sequences that violate additional, user-defined constraints.
    feasibilityVector = ones(1,M(idx));
    if strcmp(additionalConstraints, 'no')
    else
        for mm = 1:M(idx)
            feasibilityVector(mm) = feasibilityFunction(availableSequencesTemp(:,mm), uOpt, idx);
            % feasibilityFunctionIDSFID must be defined by the
            % user. It classifies feasible (function Output 1) and
            % infeasible (function Output 0) new sequences. uOpt - so far
            % optimized signal; idx current dimension that is optimized;
            % availableSequencesTemp(:,mm) - sequence checked for
            % feasibility
        end
    end

    % Remove sequences that are visited more than others and sequences that
    % are infeasible
    availableSequences = availableSequencesTemp(:, (visitsLevels{idx}==min(visitsLevels{idx}))' & feasibilityVector==1);

    % Violate the Latin hypercube if there are no sequences that are feasible
    % AND fulfill the Latin hypercube criterium.
    if isempty(availableSequences)
        messageWarning = ['To continue the optimization for input number ', num2str(idx), ' the latin hypercube criterion is violated.'];
        warning(messageWarning);
        kk = min(visitsLevels{idx});
        % Employ the sequence that is feasible with the fewest previous
        % visits.
        while isempty(availableSequences)
            availableSequences = availableSequencesTemp(:, (visitsLevels{idx}==kk)' & feasibilityVector==1);
            kk = kk+1;
            if kk == max(visitsLevels{idx})+1
                break
            end
        end
    end

    % Stop the function, if there are no feasible sequences at all.
    if isempty(availableSequences)
        error('No feasible sequences. Consider redefining the additional constraints.')
    end

    % Adjust sequences that exceed the velocity constraints.
    for ll = 1 : size(availableSequences, 2)
        kk=1;
        if availableSequences(1,ll)>uOpt{idx}(end)
            while (uOpt{idx}(end) + kk*velocityConstraints(2,idx)<availableSequences(kk,ll))
                availableSequences(kk,ll) = uOpt{idx}(end) + velocityConstraints(2,idx)*kk;
                kk=kk+1;
                if kk > Lmax(idx)
                    break
                end
            end
        elseif availableSequences(1,ll)<uOpt{idx}(end)
            while uOpt{idx}(end) + kk*velocityConstraints(1,idx)>availableSequences(kk,ll)
                availableSequences(kk,ll) = uOpt{idx}(end) + velocityConstraints(1,idx)*kk;
                kk=kk+1;
                if kk > Lmax(idx)
                    break
                end
            end
        end
    end
    % Transform possible new sequences to proxy regressor space.
    tIteration = (Ts:Ts:(numberOptDataPoints(idx)+Lmax(idx))*Ts)';
    uIteration = zeros(numberOptDataPoints(idx)+Lmax(idx),dim);
    for jj = 1:dim
        if jj == idx                                                                  % include already optimized data points of u(idx)
            uIteration(1:numberOptDataPoints(jj),jj) = uOpt{1,jj}(1:numberOptDataPoints(jj));
        elseif jj ~= idx
            if numberOptDataPoints(idx)+Lmax(idx) > size(cell2mat(uOpt(jj)),1)        % fill u(j) with last value if possible sequence length of u(idx) is bigger than the already optimized points of u(j)
                uIteration(1:numberOptDataPoints(jj), jj) = uOpt{1,jj}(1:numberOptDataPoints(jj));
                uIteration(numberOptDataPoints(jj)+1:numberOptDataPoints(idx)+Lmax(idx), jj) = uOpt{1,jj}(numberOptDataPoints(jj));
            elseif numberOptDataPoints(idx)+Lmax(idx) < size(cell2mat(uOpt(jj)),1)    % truncate u(j) if possible sequence length of u(idx) is smaller than the already optimied points of u(j)
                uIteration(1:numberOptDataPoints(idx)+Lmax(idx), jj) = uOpt{1,jj}(1:numberOptDataPoints(idx)+Lmax(idx));
            end
        end
    end
    XTemp = cell(1, size(availableSequences, 2));
    for ll = 1:size(availableSequences, 2)
        uIteration(numberOptDataPoints(idx)+1 : numberOptDataPoints(idx)+Lmax(idx),idx) = availableSequences(:,ll);
        yProxyTemp = lsim(Sys,uIteration,tIteration);
        XTemp{ll} = [uIteration, yProxyTemp];
    end
    % Calculate J1 for each available sequence and for all available lengths in proxy regressor space.
    [J1, minDist] = deal(zeros(Lmax(idx), size(availableSequences, 2)));
    XOptIteration = XTemp{1}(1:numberOptDataPoints(idx),:);                       % already optimized data points in this iteration.
    for ll = 1:size(availableSequences, 2)
        minDist(:,ll) = min(pdist2(XOptIteration, XTemp{ll}(numberOptDataPoints(idx)+1:numberOptDataPoints(idx)+Lmax(idx), :)), [], 1);
        J1(:,ll) = cumsum(minDist(:,ll));
    end
    % Calculate Factor F for each available sequence.
    F = zeros(1, size(availableSequences, 2));
    for ll = 1 : size(availableSequences, 2)
        XnewTemp = XTemp{ll}(numberOptDataPoints(idx)+1:numberOptDataPoints(idx)+Lmax(idx),:);
        innerDistances = zeros(Lmax(idx)-1, 1);
        for kk = 1:Lmax(idx)-1
            innerDistances(kk) = pdist2(XnewTemp(kk,:), XnewTemp(kk+1,:));
        end
        F(ll) = mean(innerDistances);
    end
    F = F.^-1;
    % Calculate the quality function from J1, J2, and F.
    qualityFunction = J1 - lambda(idx) * F.*repmat(J2{idx}, 1, size(availableSequences,2));
    % Get the optimal amplitudelevel and the corresponding optimal sequence
    % length.
    [maximumQualityFunction, optSeqLengths] = max(qualityFunction);
    [~, idxAmplitudeLevel] = max(maximumQualityFunction);
    idxSeqLength = optSeqLengths(idxAmplitudeLevel);
    % Concatenate sequence to existing signal.
    uOpt{idx} = [uOpt{idx}; availableSequences(1:idxSeqLength, idxAmplitudeLevel)];
    % Get the nearest neighbor amplitude level of the last u value to
    % determine the 'real' amplitude level that was selected.
    lastU = uOpt{idx}(numberOptDataPoints(idx)+idxSeqLength);
    idxRealAmplitudeLevel = knnsearch(levels{idx},lastU);
    visitsLevels{idx}(idxRealAmplitudeLevel) = visitsLevels{idx}(idxRealAmplitudeLevel)+1;
    iterations(1,idx) = iterations(1,idx)+1;
    % Remove the data point [0, 0] that is only used to start the optimization.
    if all(existingData(1,:)== 0) && size(existingData,1) == 1 && mean(iterations) == 2
        for jj=1:dim
            uOpt{jj}(1) = [];
        end
    end
end

%% Data Adjustments
[uProxy, u] = deal(zeros(N, dim));
for jj=1:dim
    uProxy(1:N,jj) = uOpt{jj}(1:N);
    u(:,jj) = rescale(uProxy(:,jj), amplitudeConstraints(1,jj), amplitudeConstraints(2,jj), "InputMin", 0, "InputMax", 1);
end
t = Ts : Ts : N*Ts;
yProxy = lsim(Sys,uProxy,t);

for jj = 1:dim
    if isempty(find(visitsLevels{jj}==0, 1))
    else
        messageWarning = ['Input number ', num2str(jj), ' exhibits an underpopulated Latin hypercube. Consider reducing the number of visitable amplitude levels M.'];
        warning(messageWarning);
    end
end

end