function feasibility = feasibilityFunction(currentSequence,uOpt,idx)
% This code demonstrates an instance where additional constraints are 
% incorporated through the feasibility function.

%% Construct current sequence in combination with already optimized data.
dim = size(uOpt, 2);
numberOptDataPoints = cell2mat(cellfun(@length, uOpt, 'UniformOutput', false));
Lmax = length(currentSequence);


 u = zeros(Lmax,dim);
    for jj = 1:dim
        if jj == idx                                                                  % include already optimized data points of u(idx)
            u(:,jj) = currentSequence;
        elseif jj ~= idx
            if numberOptDataPoints(idx)+Lmax > size(cell2mat(uOpt(jj)),1)        % fill u(j) with last value if possible sequence length of u(idx) is bigger than the already optimized points of u(j)
                u(:,jj) = [uOpt{1,jj}(numberOptDataPoints(idx)+1 : numberOptDataPoints(jj)); ones(numberOptDataPoints(idx)+Lmax-numberOptDataPoints(jj),1)*uOpt{1,jj}(numberOptDataPoints(jj))]; 
            elseif numberOptDataPoints(idx)+Lmax < size(cell2mat(uOpt(jj)),1)    % truncate u(j) if possible sequence length of u(idx) is smaller than the already optimied points of u(j)
                u(:,jj) = uOpt{1,jj}(numberOptDataPoints(idx)+1:numberOptDataPoints(idx)+Lmax);
            end
        end
    end

%% Check for feasibility.
feas = zeros(Lmax,1);
for ii = 1 : Lmax
    feas(ii) = u(ii,2)<=u(ii,1); % input 2 needs to be smaller than input 1 in normalized proxy regressor space.
end

feasibility = all(feas);

end