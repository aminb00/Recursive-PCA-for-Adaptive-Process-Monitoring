function s = sturm_count(alpha, beta, l)
% STURM_COUNT  Counts eigenvalues of Gamma strictly less than l
%
% INPUT:
%   alpha - diagonal of Gamma (m1 x 1)
%   beta  - subdiagonal of Gamma (m1-1 x 1)
%   l     - evaluation point
%
% OUTPUT:
%   s - number of eigenvalues of Gamma strictly less than l

    m1 = length(alpha);
    
    % Initialize the sequence
    % p(0) = 1
    % p(1) = alpha(1) - l
    % p(v) = (alpha(v) - l)*p(v-1) - beta(v-1)^2 * p(v-2)
    
    p_prev = 1;               % p_{v-2}  (p_0)
    p_curr = alpha(1) - l;   % p_{v-1}  (p_1)
    
    % Count sign changes
    s = 0;
    
    % Sign change between p_0 and p_1
    % Convention: if p_curr = 0, sign = -sign(p_prev)
    if sign_sturm(p_curr) ~= sign_sturm(p_prev)
        s = s + 1;
    end
    
    % --- Recurrence for v = 2, ..., m1 ---
    for v = 2:m1
        
        p_next = (alpha(v) - l) * p_curr - beta(v-1)^2 * p_prev;
        
        % Sign change between p_{v-1} and p_v
        if sign_sturm(p_next) ~= sign_sturm(p_curr)
            s = s + 1;
        end
        
        % Update for next step
        p_prev = p_curr;
        p_curr = p_next;
        
    end

end

% Helper function for sign 
% Sturm convention: if p = 0, use opposite sign of predecessor
% We implement as: sign(0) = -1 (standard Sturm convention)
function sg = sign_sturm(p)
    if p > 0
        sg = 1;
    elseif p < 0
        sg = -1;
    else
        sg = -1;  % convention: zero counts as negative
    end
end