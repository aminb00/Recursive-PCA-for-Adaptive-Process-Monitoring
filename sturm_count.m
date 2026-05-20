function s = sturm_count(alpha, beta, l)
% STURM_COUNT  Conta gli autovalori di Gamma strettamente minori di l
% usando la Sturm sequence (Li et al. 2000, Eq. 32)
%
% INPUT:
%   alpha - diagonale di Gamma (m1 x 1)
%   beta  - sottodiagonale di Gamma (m1-1 x 1)
%   l     - punto di valutazione
%
% OUTPUT:
%   s - numero di autovalori di Gamma strettamente minori di l

    m1 = length(alpha);
    
    % --- Inizializza la sequenza ---
    % p(0) = 1
    % p(1) = alpha(1) - l
    % p(v) = (alpha(v) - l)*p(v-1) - beta(v-1)^2 * p(v-2)
    
    p_prev = 1;               % p_{v-2}  (p_0)
    p_curr = alpha(1) - l;   % p_{v-1}  (p_1)
    
    % Conta cambi di segno
    s = 0;
    
    % Cambio di segno tra p_0 e p_1
    % Convenzione: se p_curr = 0, segno = -segno(p_prev)
    if sign_sturm(p_curr) ~= sign_sturm(p_prev)
        s = s + 1;
    end
    
    % --- Ricorrenza per v = 2, ..., m1 ---
    for v = 2:m1
        
        p_next = (alpha(v) - l) * p_curr - beta(v-1)^2 * p_prev;
        
        % Cambio di segno tra p_{v-1} e p_v
        if sign_sturm(p_next) ~= sign_sturm(p_curr)
            s = s + 1;
        end
        
        % Aggiorna per il passo successivo
        p_prev = p_curr;
        p_curr = p_next;
        
    end

end

% --- Funzione helper per il segno ---
% Convenzione Sturm: se p = 0, usa segno opposto al predecessore
% Implementiamo come: sign(0) = -1 (convenzione standard per Sturm)
function sg = sign_sturm(p)
    if p > 0
        sg = 1;
    elseif p < 0
        sg = -1;
    else
        sg = -1;  % convenzione: zero conta come negativo
    end
end