%Minimal implementation of a mehrotra lp solver
function [x,y,s,info] = mehrotra_lp_solver(A,b,c)

    warning('off','MATLAB:nearlySingularMatrix');

    fprintf('Minimal mehrotra type solver \n');
    %-----------------------------------------------------------------
    % Minimal Mehrotra LP predictor corrector
    %------------------------------------------------------------------
    [m,n] = size(A); 
    nu    = n;
    fprintf('Problem size %i constraints %i variables \n',m,n);
    %-----------------------------------
    %Initialization strategy from CVXOPT
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[zeros(n,1);b];
    x   = sol(1:n);
    
    K2  = [[speye(n), A'];[A,sparse(m,m)]];
    sol = K2\[c;zeros(m,1)];
    y   = sol(n+1:n+m);
    s   = sol(1:n);
    
    %%XXX:RANDOM INIT
    %x   = ones(n,1);
    %s   = ones(n,1);
    %y   = zeros(m,1);
    
    tau   = 1;
    kappa =1;
    clear sol
    
    if(min(x)<0) %if x is not feasible shift it into feasibility
        a = min(x);
        x = (1-a)*ones(n,1)+x;
    end
    
    if(min(s)<0) %if s is not feasible shift it into feasibility
        a = min(s);
        s = (1-a)*ones(n,1)+s;
    end
    
    %------Calculate the initial mu
    mu = (s'*x+tau*kappa)/(nu+1);
    mu0 = mu;
    
    %--------------------
    %Evaluate the residuals
    rp = tau*b-A*x;
    rd = -tau*c+A'*y+s;
    rg = kappa-b'*y+c'*x;
    
    %Evaluate the residual norms
    nrp = norm(rp);
    nrd = norm(rd);
    nrg = norm(rg);
    
    % Save the original norms
    nrp0 = nrp;
    nrd0 = nrd;
    nrg0 = nrg;
    
    gap  = c'*x-b'*y;
    fprintf('%2i a       %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',0,nan,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,nan);
    
    %Algorithm parameters
    max_iter = 100;
    
    %--------------------------------------
    %Main SOCP solver iteration
    for iter=1:max_iter
       
        %Evaluate the hessian and scaled variables
        %Evaluate the scaling point 
        H = diag(sparse(s./x));
    
        %Build the rhs term 
        %Solve the affine scalin g direction
        [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-s,-kappa,[]);
    
        ratios = [-x./dx;-s./ds;-tau/dtau;-kappa/dkappa;1];
        a_max  = min(ratios(find(ratios>0)));
    
        %Now calculate the centering direction
        %-------------------------------------
        sigma = (1-a_max)^3;
    
        %Build the second order correction term, this should have a better form
        so       = (1-sigma)*dx.*ds./x;
        kt_so    = (1-sigma)*dtau*dkappa/tau;
          
        %so    = 0;
        %kt_so = 0;
        
        %Solve the affine scaling direction
        [dy_c,dx_c,dtau_c,ds_c,dkappa_c,residual_norm_c,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,...
                                                                                                (1-sigma)*rp,...
                                                                                                (1-sigma)*rd,...
                                                                                                (1-sigma)*rg,...
                                                                                                -s    +sigma*mu*1./x  -so,...
                                                                                                -kappa+sigma*mu*1/tau -kt_so,...
                                                                                                slv_aug);
    
    
        %Decrement 
        ratios = [-x./dx_c;-s./ds_c;-tau/dtau_c;-kappa/dkappa_c;1];
        a_max  = min(ratios(find(ratios>0)));
        a       = a_max*0.98; 
        % a       = min(a,1); 
        
        x_old = x;
        y_old = y;
        s_old = s;
        t_old = tau;
        k_old = kappa;
    
        y         = y+a*dy_c;
        x         = x+a*dx_c;
        s         = s+a*ds_c;
        tau       = tau+a*dtau_c;
        kappa     = kappa+a*dkappa_c;
    
        mu        = x'*s+tau*kappa;
        mu        = mu/(nu+1);
        
        %Feasiblity check
        x_slack = min(x);
        s_slack = min(s);
        if(x_slack<0||s_slack<0||tau<0||kappa<0)
            fprintf('Infeasible centering step \n');
            return;
        end
    
        %--------------------
        %Evaluate the residuals
        rp = tau*b-A*x;
        rd = -tau*c+A'*y+s;
        rg = kappa-b'*y+c'*x;
    
        nrp   = norm(rp);
        nrd   = norm(rd);
        nrg   = norm(rg);
        
        gap   = (mu*(nu+1)-tau*kappa)/tau;
        
        fprintf('%2i a %3.3e s %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',...
                iter,a,sigma,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,residual_norm_c);

        if(nrp/nrp0 < 1.e-8 && nrd/nrd0 < 1.e-8 && mu/mu0 < 1.e-8)
            break;
        end
    end %end main loop

    x = x/tau;
    s = s/tau;
    y = y/tau;
    info = struct;
    info.iter = iter;
end

