function out = numSources(D)
    % D is the diagonal eigen value matrix of the matrix of empirical
    % covariance of the antenna data
    num = NaN;
    for i = size(D, 1):-1:3
        
        if (D(i-1, i-1)-D(i-2, i-2))/(D(i, i)-D(i-1, i-1))<=0.1
            num = i-1;
        end
        
    end
    
    out = num;
end