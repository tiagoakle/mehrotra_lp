%Minimal implementation of a mehrotra lp solver

% the opts structure contains the following options 
% opts.max_iter % maximum number of iterations 
% opts.centrality one of none functional 2norm 
%                 none uses a long-step method 
%                 functional keeps nu log (x's) + f(x) + f(s) - nulog(nu) + nu < centrality
%                 2norm keeps ||Xz-mu||_2 \leq centrality
% opts.ini_mehrotra true uses merhotras initialization (defaults to false if cetnrality is not none)

function [x,y,s,info] = mehrotra_lp_solver(A,b,c,opts) 
   
    warning('off','MATLAB:nearlySingularMatrix');

    
    if(opts.verbose) fprintf('Minimal mehrotra type solver \n'); end;
    %-----------------------------------------------------------------
    % Minimal Mehrotra LP predictor corrector
    %------------------------------------------------------------------
    [m,n] = size(A); 
    nu    = n;
    if(opts.verbose) fprintf('Problem size %i constraints %i variables \n',m,n); end;
    %-----------------------------------
    
    if(nargin == 3)
        opts = struct;
    end
    opts = get_default_options(opts);
    
    if(~opts.use_centrality)
        bk_iter_predictor = 0;
        bk_iter_corrector = 0;
    end 
    if(opts.verbose) fprintf('\t Centrality %s, \n \t Initialization Mehrotra %g \n',opts.centrality_type,opts.ini_mehrotra); end 
    if(opts.ini_mehrotra)
        %Initialization strategy from CVXOPT
        K2  = [[speye(n), A'];[A,sparse(m,m)]];
        sol = K2\[zeros(n,1);b];
        x   = sol(1:n);
        
        K2  = [[speye(n), A'];[A,sparse(m,m)]];
        sol = K2\[c;zeros(m,1)];
        y   = sol(n+1:n+m);
        s   = sol(1:n);
    else
        x = ones(n,1);
        s = ones(n,1);
        y = zeros(m,1);
    end
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
    if(opts.verbose)
        fprintf('%2i a       %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e\n',0,nan,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,nan);
    end 
    
    %In case we use centrality backtrack constant 
    bk_constant = 0.8;
    max_bk_iter     = 100;

    
    %--------------------------------------
    %Main SOCP solver iteration
    for iter=1:opts.max_iter
       
        %Evaluate the hessian and scaled variables
        %Evaluate the scaling point 
        H = diag(sparse(s./x));
    
        %Build the rhs term 
        %Solve the affine scalin g direction
        [dy,dx,dtau,ds,dkappa,res_norm,slv_aug] = linear_solver(H,tau,kappa,A,b,c,m,n,rp,rd,rg,-s,-kappa,[]);
    
        ratios = [-x./dx;-s./ds;-tau/dtau;-kappa/dkappa;1];
        a_max  = min(ratios(find(ratios>0)));
        
        bk_iter_predictor = 0;
        if(opts.use_centrality) %Find the largest alpha that keeps ||Xz-mu||^2_2 < centality
            for bk_iter_predictor = 1:max_bk_iter
                cent = 0;
                if(opts.use_functional)
                    cent = potential_centrality(x+a_max*dx,s+a_max*ds,tau+a_max*dtau,kappa+a_max*dkappa);
                else
                    cent = two_norm_centrality(x+a_max*dx,s+a_max*ds,tau+a_max*dtau,kappa+a_max*dkappa,mu);
                end
                if(cent < opts.centrality)
                    break;
                else
                    a_max = a_max*bk_constant;
                end
            end
        end
    
        %Now calculate the centering direction
        %-------------------------------------
        sigma = (1-a_max)^3;
    
        %Build the second order correction term, this should have a better form
        so       = (1-sigma)*dx.*ds./x;
        kt_so    = (1-sigma)*dtau*dkappa/tau;
        
        %If no second order correction should be used then
        if(~opts.secord)
            so    = 0;
            kt_so = 0;
        end
        
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
        
        
        if(opts.use_centrality) %Find the largest alpha that keeps ||Xz-mu||^2_2 < centality
            for bk_iter_corrector = 1:max_bk_iter
                cent = 0;
                if(opts.use_functional)
                    cent = potential_centrality(x+a_max*dx_c,s+a_max*ds_c,tau+a_max*dtau_c,kappa+a_max*dkappa_c);
                else
                    cent = two_norm_centrality(x+a_max*dx_c,s+a_max*ds_c,tau+a_max*dtau_c,kappa+a_max*dkappa_c,mu);
                end
                if(cent < opts.centrality)
                    break;
                else
                    a_max = a_max*bk_constant;
                end
            end
        end

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
        if(opts.verbose)      
            fprintf('%2i a %3.3e s %3.3e pr %3.3e dr %3.3e gr %3.3e mu %3.3e gap %3.3e k/t %3.3e res_cent %3.3e bk pre %2i bk cor %2i\n',...
                iter,a,sigma,nrp/(nrp0*tau),nrd/(nrd0*tau),nrg,mu,gap,kappa/tau,residual_norm_c,bk_iter_predictor,bk_iter_corrector);
        end
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

function potential = potential_centrality(x,s,t,k)
    nu = size(x,1)+1;
    potential = nu*log(x'*s+t*k)-sum(log(x))-sum(log(s))-log(t)-log(k) - nu*log(nu);
end

function cent = two_norm_centrality(x,s,t,k,mu)
                act_res = ([x;t]).*([s;k])-mu; 
                cent = norm(act_res);
end

function opts = get_default_options(opts)
    if(~isfield(opts,'max_iter'))        opts.max_iter        = 100   ; end
    if(~isfield(opts,'centrality_type')) opts.centrality_type = 'none'; end
    if(~isfield(opts,'ini_mehrotra'))    opts.ini_mehrotra    = true  ; end
    if(~isfield(opts,'secord'))          opts.secord          = true  ; end
    if(~isfield(opts,'verbose'))         opts.verbose         = true  ; end
    
    %If we are using centrality then disable ini_mehrotra
    if(strcmp(opts.centrality_type,'none')~= 1)
        opts.ini_mehrotra = false; 
        opts.use_centrality = true;

        %if the centrality threshold is not set, then set it
        %Unset the second order correction
        if(strcmp(opts.centrality_type,'functional')==1)
            if(~isfield(opts,'centrality')) opts.centrality = 1-log(2); end
            opts.use_functional = true;
            opts.secord = false;
        elseif(strcmp(opts.centrality_type,'2norm')==1)
            if(~isfield(opts,'centrality')) opts.centrality = 0.5; end
            opts.use_functional = false;
        end
    else
        opts.use_centrality = false;
    end

end
