A = rand(100);

maxerror = 0;

for i = 2:100
    
    A = rand(6);

    [L, U, P] = LUPDecomposition(A);

    error = norm(P*A-L*U);
    
    if error > maxerror
        maxerror = error;
    end
end

maxerror
    