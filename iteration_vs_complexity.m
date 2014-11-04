%Create a sequence of lp problems of increasing size 
%solve each problem with a predictor-corrector with short-step 
%and with a predictor-corrector with no centrality constraints

sizes = [100:100:200];
samples = 3;

results = zeros(length(sizes)*samples,4);
result_ix = 1;
for size_ix = 1:length(sizes)
    for sample = 1:samples
        n = sizes(size_ix);
        m = ceil(n/10);
        A = sprandn(m,n,0.8);
        b = A*rand(n,1);
        c = A'*randn(m,1) + rand(n,1);

        %No centrality and ones for initialization
        opts = struct;
        opts.centrality_type = 'none';
        opts.ini_mehrotra    = false;
        opts.secord          = false;
        [x,y,s,info_none] = mehrotra_lp_solver(A,b,c,opts);
      
        %No centrality and ones for initialization
        opts = struct;
        opts.centrality_type = 'functional';
        opts.secord          = false;
        [x,y,s,info_func] = mehrotra_lp_solver(A,b,c,opts);

        %Centrality and ones for initialization
        opts = struct;
        opts.centrality_type = '2norm'
        opts.secord          = false;
        [x,y,s,info_cent] = mehrotra_lp_solver(A,b,c,opts);
        results(result_ix,:) = [n, info_none.iter, info_func.iter, info_cent.iter];
        result_ix = result_ix + 1;
    end
end

%Average over the samples and create a new results 
avg_results = zeros(length(sizes),4);
for size_ix = 1:length(sizes)
    n = sizes(size_ix);
    ixs = find(results(:,1)==n);
    avg_row = mean(results(ixs,2:end));
    avg_results(size_ix,:) = [n avg_row];
end

csvwrite('AverageIterationVsComplexity.csv',avg_results);


