function lambda_i = bisection_eig(alpha, beta, i, eps)
% BISECTION_EIG  Trova l'i-esimo autovalore più grande di Gamma
% tramite bisection + Sturm sequence
% Following Li et al. (2000), Section 4
%
% INPUT:
%   alpha    - diagonale di Gamma (m1 x 1)
%   beta     - sottodiagonale di Gamma (m1-1 x 1)
%   i        - indice autovalore (1 = più grande, 2 = secondo, ...)
%   eps      - precisione (default 1e-10)
%
% OUTPUT:
%   lambda_i - i-esimo autovalore più grande di Gamma

    if nargin < 4
        eps = 1e-10;
    end

    m1 = length(alpha);

    % --- Bounds di Gershgorin ---
    % Ogni autovalore di Gamma sta in [l_min, l_max]
    % Per una matrice tridiagonale simmetrica:
    % riga 1:    [alpha(1) - |beta(1)|,          alpha(1) + |beta(1)|]
    % riga j:    [alpha(j) - |beta(j-1)| - |beta(j)|, alpha(j) + ...]
    % riga m1:   [alpha(m1) - |beta(m1-1)|,      alpha(m1) + |beta(m1-1)|]

    % Calcola il raggio di Gershgorin per ogni riga
    gershgorin_radius = zeros(m1, 1);
    gershgorin_radius(1) = abs(beta(1));
    for j = 2:m1-1
        gershgorin_radius(j) = abs(beta(j-1)) + abs(beta(j));
    end
    gershgorin_radius(m1) = abs(beta(m1-1));

    l_max = max(alpha + gershgorin_radius) + eps;
    l_min = min(alpha - gershgorin_radius) - eps;

    % --- Bisection ---
    % Vogliamo lambda_i tale che:
    % s(lambda_i) = m1 - i
    % cioè: esattamente (m1 - i) autovalori sono strettamente minori di lambda_i
    % (perché ci sono i-1 autovalori più grandi di lambda_i)

    target = m1 - i;  % numero di autovalori che devono stare sotto lambda_i

    a = l_min;
    b = l_max;

    max_iter = ceil(log2((b - a) / eps)) + 10;

    for iter = 1:max_iter

        l_mid = (a + b) / 2;

        s_mid = sturm_count(alpha, beta, l_mid);

        if s_mid <= target
            % lambda_i è a destra di l_mid (o è l_mid stesso)
            a = l_mid;
        else
            % lambda_i è a sinistra di l_mid
            b = l_mid;
        end

        % Criterio di convergenza
        if (b - a) < eps
            break
        end

    end

    lambda_i = (a + b) / 2;

end