function [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,r1,r2,r3,r5,r4,factorization)
    
    %regularizations 
    delta = 1.e-10;
    gamma = 1.e-10;
    max_iter_ref_rounds = 10;

    nf = 0;
    print_t = false;
    %Check if there is an existing factorization
    if(~isempty(factorization))
        slv_aug = factorization;
    else
    %No factorization provided, factor the matrix
     slv_aug           = factor();
    end
    
    %Eventually this will be implemented in C
    %Solves  [    A -b      ] dy   r1
    %        [-A'    c -I   ] dx   r2
    %        [b' -c'     -1 ] dt = r3
    %        [     H   I    ] ds   r5
    %        [       k    t ] dk   r4

    %Forms 
    %        [    A -b      ] dy   r1
    %        [-A'    c -I   ] dx   r2
    %        [b' -c'     -1 ] dt = r3
    %        [     H   I    ] ds   r5
    %        [       h   1  ] dk   r6
    %with h = k/t r6 = r4/tau

    h  = kappa/tau;
    r6 = r4; %This variable is not necessary

    %Forms 
    %        [    A  -b      ] dy     r1
    %        [A' - H -c      ] dx     r7
    %        [b' -c'  h      ] dt =   r8
    %        [    H     I    ] ds     r5
    %        [        h   1  ] dk     r6
    % with r7 = -r2-r5; r8 = r3+r6;
    
    r7 = -(r2+r5);
    r8 = r3+r6;

    %  
    %        [    A  -b     ] dy     r1
    %        [A'  -H -c     ] dx     r7
    %        [       h_2    ] dt =   r9
    %        [    H     I   ] ds     r5
    %        [       h    1 ] dk     r6
     
    %We let 
    %[b' -c'] [     A]^{-1}[-b]     =  h_1b
    %         [A'  -H]     [-c] 

    %[b' -c'] [     A]^{-1}[r1]     =  r_7b
    %         [A'  -H]     [r7] 
    
    %Now calculate r_7b, h_1b. Solve 
    % tm_1 =  [     A]^{-1}[b ]    
    %         [A'  -H]     [-c] 
    % and keep tm_1 it to calculate h_1b and r_7b more efficiently
   
    if(print_t) fprintf('Norm [b;-c] %g\n',norm([b;-c])); end;

    tm_1   = slv_aug([b;-c]); 
  
    if(print_t) fprintf('Norm tm1 %g\n',norm(tm_1)); end;
    
    h_1b   = tm_1'*[-b;-c];
    r_7b   = tm_1'*[r1;r7];
 
    r9  = r8 - r_7b;
    h_2 = h  - h_1b;
    
    %We now start the back substitution
    dt  = r9/h_2;

    %reuse tm_1
    tm_1 = slv_aug([r1;r7]-dt*[-b;-c]);

    %Extract dy dx from dt
    dy   = tm_1(1:m);
    dx   = tm_1(m+1:m+n);

    %Now back substitute to get ds and dkappa
    ds   = -A'*dy+dt*c-r2;
    dk   = r6-h*dt;
    
    %Check the residuals
    n_res_1 = norm(A*dx-dt*b-r1);
    n_res_2 = norm(-A'*dy + dt*c - [zeros(nf,1);ds] - r2);
    n_res_3 = norm(b'*dy-c'*dx -dk-r3);
    n_res_5 = norm(H*dx+ds-r5);
    n_res_4 = norm(kappa*dt+tau*dk-r4);
    
    res_norm = norm([n_res_1,n_res_2,n_res_3,n_res_4,n_res_5]);

    if(print_t) fprintf('Residuals of mixed C solve r1 %g, r2 %g, r3 %g, r5 %g, r4 %g \n',n_res_1,n_res_2,n_res_3,n_res_5,n_res_4); end;
    
    dx     = dx; 
    dtau   = dt;
    dy     = dy;
    ds     = ds;
    dkappa = dk;
    
         
    function h_solver = factor()
       %Factorizes
        %[dI  A    ]    = PLDL'P'
        %[A' -H -gI]
        [L,D,P] = ldl([[delta*speye(m,m),A];[A',-H-gamma*speye(n,n)]]);
        %Instantiate a handle to the solver
        h_solver = @(y)solve(L,D,P,y);
        if(print_t) fprintf('Sizes m: %i n %i \n',m,n); end
        linopts_lt = struct;
        linopts_ut = struct;
        linopts_lt.LT = true;
        linopts_ut.UT = true;
    
        d = diag(D);
        if(print_t) fprintf('min(d): %g, max(d): %g \n',full(min(d)),full(max(d))); end;
        if(print_t) fprintf('min diag L: %g, max |L|: %g \n', min(abs(full(diag(L)))),full(max(max(abs(L))))); end;
    
        function x = solve(L,D,P,y)
           
            %Calculate the norms to stop the iterative refinement rounds
            n_y = norm(y,'inf');
            
            %For debug 
            y_t = y;
            %PLDL'P'x = b 
            if(print_t) fprintf('Norm rhs: %g\n', norm(y)); end;
            x  = (P*(L'\(D\(L\(P'*y)))));
            
            %Iterative refinement rounds
            %
            %[      A ] x(1)   = y(1)
            %[A'   -H ] x(2)     y(2)
            for(iter_round = 1:max_iter_ref_rounds) 
                ri1 = y(1:m) - A*x(m+1:m+n);
                ri2 = y(m+1:m+n) - A'*x(1:m) + H*x(m+1:m+n);
                n_res = norm([ri1;ri2],'inf');
                if(print_t) fprintf('Round %i: res norm %g,  relative residual norm %g\n',iter_round,n_res,n_res/n_y); end
                if(n_res/n_y<1.e-15)
                    break;
                end

            %Calculate the correction
            x  = x + (P*(L'\(D\(L\(P'*[ri1;ri2])))));
            end %End of iterative refinement loop

            %Print the final residuals
            if(print_t)
                ri1 = y(1:m) - A*x(m+1:m+n);
                ri2 = y(m+1:m+n) - A'*x(1:m) + H*x(m+1:m+n);
                n_res = norm([ri1;ri2],'inf');
                fprintf('Final iterative ref  res norm %g, relative residual norm %g \n',n_res,n_res/n_y); 
            end
            
        end
    
    end
 
end



  
