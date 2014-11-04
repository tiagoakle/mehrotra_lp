
m = 10;
n = 100;
A = sprandn(m,n,0.8);
b = A*rand(n,1);
c = A'*randn(m,1) + rand(n,1);

[x,y,s,info] = mehrotra_lp_solver(A,b,c);
%No centrality and ones for initialization
[x,y,s,info] = mehrotra_lp_solver(A,b,c,0,false);
%Centrality and ones for initialization
[x,y,s,info] = mehrotra_lp_solver(A,b,c,1,false);

