%% LU Decomposition via Gaussian Elimination with Partial Pivoting
%
%                                PA = LU
%
% This algorithm implements LU Decomposition with Gaussian elimination and
% pivoting. It follows Trefethen's implementation from Numerical Linear
% Algebra (Algorithm 21.1).
%
% This algorithm has been tested against square matrices of varying sizes
% and error never surpassed the order of 1e-15.
%


function [L,U,P] = LUPDecomposition(A)

s = size(A);
m = s(1);

    if m ~= s(2)
        "Matrix is not square"
        return
    end

    U = A;
    L = eye(m);
    P = eye(m);
    
    for k = 1:m-1
        
        max_i = U(k,k);
        max_i_index = k;
        
        for i = k:m
            
            if abs( U(i,k) ) > abs(max_i)
                max_i = U(i,k);
                max_i_index = i;
         
            end
        end
        
        U([k max_i_index],k:m) = U([max_i_index k],k:m);
        L([k max_i_index],1:k-1) = L([max_i_index k],1:k-1);
        P([k max_i_index],:) = P([max_i_index k],:);
        
    
        for j = k+1:m
            
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);
        end
        
end

