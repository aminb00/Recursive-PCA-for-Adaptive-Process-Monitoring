function lambda_i = bisection_eig(alpha, beta, i, eps)
% BISECTION_EIG  Finds the i-th largest eigenvalue of Gamma
% via bisection + Sturm sequence
% Following Li et al. (2000), Section 4
%
% INPUT:
%   alpha    - diagonal of Gamma (m1 x 1)
%   beta     - subdiagonal of Gamma (m1-1 x 1)
%   i        - eigenvalue index (1 = largest, 2 = second largest, ...)
%   eps      - precision (default 1e-10)
%
% OUTPUT:
%   lambda_i - i-th largest eigenvalue of Gamma

    if nargin < 4
        eps = 1e-10;
    end

    m1 = length(alpha);

    %Gershgorin bounds

    % Every eigenvalue of Gamma lies in [l_min, l_max]
    % For a symmetric tridiagonal matrix:
    % row 1:    [alpha(1) - |beta(1)|,          alpha(1) + |beta(1)|]
    % row j:    [alpha(j) - |beta(j-1)| - |beta(j)|, alpha(j) + ...]
    % row m1:   [alpha(m1) - |beta(m1-1)|,      alpha(m1) + |beta(m1-1)|]

    %Calculate Gershgorin radius for each row
    gershgorin_radius = zeros(m1, 1);
    gershgorin_radius(1) = abs(beta(1));
    for j = 2:m1-1
        gershgorin_radius(j) = abs(beta(j-1)) + abs(beta(j));
    end
    gershgorin_radius(m1) = abs(beta(m1-1));

    l_max = max(alpha + gershgorin_radius) + eps;
    l_min = min(alpha - gershgorin_radius) - eps;

   %Bisection 

    % We want lambda_i such that:
    % s(lambda_i) = m1 - i
    % i.e.: exactly (m1 - i) eigenvalues are strictly less than lambda_i
    % (because there are i-1 eigenvalues larger than lambda_i)

    target = m1 - i;  % number of eigenvalues that must be below lambda_i

    a = l_min;
    b = l_max;

    max_iter = ceil(log2((b - a) / eps)) + 10;

    for iter = 1:max_iter

        l_mid = (a + b) / 2;

        s_mid = sturm_count(alpha, beta, l_mid);

        if s_mid <= target
            % lambda_i is to the right of l_mid (or is l_mid itself)
            a = l_mid;
        else
            % lambda_i is to the left of l_mid
            b = l_mid;
        end

        % Convergence criterion
        if (b - a) < eps
            break
        end

    end

    lambda_i = (a + b) / 2;

end