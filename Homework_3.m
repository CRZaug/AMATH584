%% Problem 1
% Here we will develop a numerical algorithm that implements a modified
% Gram-Scmidt orthogonalization procedure to perform QR decomposition. We
% compare it to qrfactor.m (which uses Householder Triangularization)
% and MATLAB's QR decomposition algorithm.
%
% We perform QR decomposition of 50 different matrices A of size $m\times n$ where $m> n$,
% studying the impact of condition number on error.
 
cA = [];
mgs = [];
qrf = [];
qrm = [];

for i = 1:50
    A = rand(100,50);

    condA = cond(A);

    [Q,R] = modifiedGramSchmidtQR(A);
    [Q1,R1] = qrfactor(A);
    [Q2, R2] = qr(A);

    [norm(Q*R-A), norm(Q1*R1-A), norm(Q2*R2-A)];
    
    cA = [cA, condA];
    mgs = [mgs, norm(Q*R-A)];
    qrf = [qrf, norm(Q1*R1-A)];
    qrm = [qrm,norm(Q2*R2-A)];
    

end

maxConditionNumber = max(cA)

fig1 = figure(1);

subplot(3,1,1);
plot(cA, mgs,'.');
title("i)");
ylabel("Modified Gram-Schmidt Error");


subplot(3,1,2);
plot(cA, qrf,'.');
title("ii)");
ylabel("qrfactor.m Error");

subplot(3,1,3);
plot(cA, qrm,'.');
title("iii)");
ylabel("MATLAB QR Algorithm Error");

han=axes(fig1,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
xlabel(han,'cond(A)');
title({"Comparision of QR factorization algorithms";' '});


%% 
% It appears that all three algorithms do well regardless of the conditioning of the
% matrix, where the error produced by each algorithm is given by $$ E = ||\mathbf{QR}-\mathbf{A}||.$
%
% The error that each algorithm produces is quite small, though Modified
% Gram-Schmidt produces errors an order of magnitude larger than qrfactor.m
% and MATLAB's QR decomposition algorithm, at 1e-15 compared to 1e-14,
% respectively. 
%
% However, we notice that the condition numbers in the above figure
% are all relatively close to 1. Let us see how the algorithms do on an
% extremely poorly conditioned matrix. Below, we perform all three algorithms 
% 50 times on an ill-conditioned
% matrix (two columns are exact copies of one another) and average the
% resulting errors produced by each algorithm.

clear all;

cA = [];
mgs = [];
qrf = [];
qrm = [];

for i = 1:50
    A = rand(100,99);
	A = [A, A(:,1)];

    condA = cond(A);

    [Q,R] = modifiedGramSchmidtQR(A);
    [Q1,R1] = qrfactor(A);
    [Q2, R2] = qr(A);
    
    cA = [cA, condA];
    mgs = [mgs, norm(Q*R-A)];
    qrf = [qrf, norm(Q1*R1-A)];
    qrm = [qrm, norm(Q2*R2-A)];
    

end

cAmean = mean(cA) % The mean condition number
results = [mean(mgs),mean(qrf), mean(qrm)];

    
fig2 = figure(2);

bar(results);
set(gca, 'xticklabel',["Modified GS"; "qrfactor.m"; "MATLAB QR"]);
ylabel("Error")

%%
% This graph was created by creating 50 ill-conditioned matrices with an
% average condition number above 1e16. The Modified Gram-Schmidt algorithm
% had the smallest error by an order of magnitude compared to MATLAB's
% algorithm and qrfactor.m. The error is not significantly larger than the
% errors produced by the well-conditioned matrices in our first trial.
%
%% Problem 2
% We will study what occurs when plotting a tiny interval $x\in[1.920, 20.80]$ 
% of the polynomial $p(x) = (x-2)^9$ vs its expansion. We will use the step
% size $\delta x = 0.001$. 

a = 1.920;
b = 2.080;
dx = 0.001;

x = [a:dx:b];

fig3 = figure(3);
subplot(1,2,1)
plot(x,p1(x))
title("a) Expansion of p(x)")

subplot(1,2,2)
plot(x,p2(x))
title("b) p(x) (not expanded)")

%%
% Plotting the expansion of $p(x)$ produces significantly more error than
% plotting the unexpanded form. Because the values of $p(X)$ in this
% interval
% are so small, numerical roundoff becomes an issue. When each small error is
% taken to the 9th power in the case of $p(x)$, the resulting error is much
% smaller and therefore not noticable. In the case of the expansion of $p(x)$, the numerical errors are
% taken to lower powers and therefore have a stronger cumulative effect.
% 

%% Problem 3
% In this problem, we will study the effects of the conditioning of a
% matrix under different constructions and perterbations.

%% 3a)
% Here we construct a random matrix A of size $m\times n$ where $m > n$
% and study the condition number as the size of the matrix increases.

clear all;

heatmapCA = zeros(100,100);


for m = 1:100
    newrow = zeros(1,100);
    for n = 1:100
        if m > n
            A = rand(m,n);
            condA = log(cond(A));
            newrow(n) = condA;
        end
    end
    
    heatmapCA(m,:) = newrow;
end

fig4 = figure(4);

h = heatmap(heatmapCA, 'GridVisible','off');
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar

title("Heatmap of log(cond(A)) for matrix size mxn");
xlabel("n");
ylabel("m");


%%
% Matrices with a similar number of rows and columns have higher condition
% numbers, while matrices with many fewer columns than rows have low
% condition numbers. None of the condition numbers appear unreasonably
% large, as they are all much lower than 1e16. Generally, matrices with
% more elements have larger condition numbers.
%
%% 3b)
% Now we create an ill-conditioned matrix by creating a matrix with two
% identical columns. Now not all columns are linearly independent.

A = rand(128,127);
An = [A A(:,1)];

conditionNumbers = [cond(A), cond(An)]
detAn = det(An)

%%
% The condition number of the matrix with two duplicate columns is much
% greater than the original matrix with no copies. Duplicating the column took a well-conditioned matrix
% and made it ill-conditioned. The determinant of the
% matrix with two duplicate columns is very large even though we expect it to be near 0, so the determinant is
% not a great indicator for the conditioning of our problem. This is because the
% determinant computation is unstable. The condition number is a better
% indicator of the conditioning.

%% 3c)
% Now we will add some noise to a duplicated column in a matrix and study
% what happens to the condition number. Here, we generate the noise first
% and simply increase its magnitude for $\epsilon \in [0.00001,0.01]$.
clear all;

clf;

cA = [];

A = rand(128,127);
noise = rand(128,1);
for eps = 0.00001:0.00001:.01
    
    An = [A A(:,1)+eps*noise];
    
    condA = cond(An);
    
    cA =[cA, condA];
end

fig5 = figure(5);
subplot(2,1,1);
plot(0.00001:0.00001:.01,cA,'.')
xlabel("epsilon");
ylabel("cond(A)");

subplot(2,1,2);
plot(0.00001:0.00001:.01,log(cA),'.')
xlabel("epsilon");
ylabel("log(cond(A))");

%%
% The first plot shows the condition number of the matrix as a function of
% the size of the noise, and the second shows the log of the condition
% number. We see that as the size of the noise grows, the condition number
% drops more rapidly than exponentially. Adding even small value of noise
% can greatly improve conditioning of the matrix.
%
% We will now repeat the above process, this time generating the noise
% inside the loop (so each matrix sees different noises of different
% magnitudes) and determine what happens to the condition number.
clear all;

clf;

cA = [];

A = rand(128,127);

for eps = 0.00001:0.00001:.01
    
    noise = rand(128,1);
    An = [A A(:,1)+eps*noise];
    
    condA = cond(An);
    
    cA =[cA, condA];
end

fig6 = figure(6);
subplot(2,1,1);
plot(0.00001:0.00001:.01,cA,'.');
xlabel("epsilon");
ylabel("cond(A)");

subplot(2,1,2);
plot(0.00001:0.00001:.01,log(cA),'.')
xlabel("epsilon");
ylabel("log(cond(A))");

%%
% With different random noise each time, the condition number still
% decreases faster than exponentially, but there is much noise in this
% process as well. Still, adding even small noise greatly decreases the
% condition number of a matrix and leads to a more well-conditioned
% problem.