function OUT = prewhiten(X)
    % define N and P.
    [N,P] = size(X);
    assert(N >= P);

    %SVD of covariance of X. We could also use svd(X) to proceed but N
    % can be large and so we sacrifice some accuracy for speed.
    [U,signal] = svd(cov(X));
    signal     = diag(signal);
    signal     = signal(:)';

    %Understanding which values of signal are non-zero.
    tol = eps(class(X));
    idx = (signal > max(signal)*tol);
    assert(~all(idx == 0));

    % Get the non-zero elements of our signal and corresponding columns of U.
    signal = signal(idx);
    U   = U(:,idx);

    %Compute prewhitened data.
    mu = mean(X,1);
    OUT = bsxfun(@minus,X,mu);
    OUT = bsxfun(@times,OUT*U,1./sqrt(signal));
end
