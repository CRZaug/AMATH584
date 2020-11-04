%% Modified Gram-Schmidt Algorithm

% Here we implement a modified Gram-Schmidt algorithm to perform QR
% decomposition.

% See below for implementation and analysis of this algorithm

function [Q, R] = modifiedGramScmidtQR(A)

    [m,n] = size(A);
    V = A;
    Q = zeros(m,n);
    R = zeros(n,n);
    
    for j = 1:n
        
        R(j,j) = norm(V(:,j));
        Q(:,j) = V(:,j)/norm(V(:,j));
        for k = j+1:n
            R(j,k)= (Q(:,j)'*V(:,k));
            V(:,k) = V(:,k) - R(j,k)*Q(:,j);
        end     
    end
    
end
