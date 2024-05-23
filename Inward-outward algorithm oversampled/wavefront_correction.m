function [psi_opt,phi_in,phi_out] = wavefront_correction(x,y,x_zone,y_zone,r_update,k,Z,FOM)

    [X,Y] = meshgrid(x,y); X = single(X); Y = single(Y);
    x_zone_min = min(x_zone,[],'all'); y_zone_min = min(y_zone,[],'all');
    x_zone_max = max(x_zone,[],'all'); y_zone_max = max(y_zone,[],'all');
    
    Nx = length(x);
    kx = k(:,1); ky = k(:,2); N = length(k);
    
    n_order = width(Z);
    c_in = zeros(n_order,1); c_out = c_in;
    
    % Build pre-optimization image
    psi_init = single(zeros(Nx,Nx));
    parfor ii = 1:Nx^2
        k_in_x = kx(ii); k_in_y = ky(ii);
        % Convert outputs to spatial basis
        r_update_ii = reshape(r_update(:,ii),Nx,Nx);
        r_r_k_ii = ifft(ifftshift(ifft(ifftshift(r_update_ii,1),Nx,1),2),Nx,2);
        % Convert the inputs to spatial basis
        psi = reshape(r_r_k_ii,Nx,Nx).*exp(-1i*k_in_x*X-1i*k_in_y*Y);

        psi_init = psi_init+psi;
    end
    I_init = abs(psi_init).^2;
    M = fom(I_init,FOM);
   
    I_init = reshape(I_init,Nx,Nx);
    
    figure(3)
    subplot 121
    imagesc(fliplr(interp2(I_init,3,'spline')))
    axis image
    colormap('hot')
    set(gca,'Visible','off')

    delta_M = 1;
    
    r_in = r_update; r_out = r_update.'; r_in_update = r_in; r_out_update = r_out;
    
    while delta_M > 0.1
        
        M_prev = M;

        % Build Psi_in matrix
        fprintf("Building Psi_in.\n")
        Psi_in = single(zeros(Nx^2,Nx^2));
        parfor ii = 1:Nx^2
            k_in_x = kx(ii); k_in_y = ky(ii);
            % Convert outputs to spatial basis
            r_update_ii = reshape(r_in_update(:,ii),Nx,Nx);
            r_r_k_ii = ifft(ifftshift(ifft(ifftshift(r_update_ii,1),Nx,1),2),Nx,2);
            % Convert the inputs to spatial basis
            psi = reshape(r_r_k_ii,Nx,Nx).*exp(-1i*k_in_x*X-1i*k_in_y*Y);
            
            Psi_in(:,ii) = psi(:);
        end
        
        
        % Optimize inputs
        fprintf("Optimizing inputs.\n")
        opt.algorithm = NLOPT_LD_LBFGS; % Choose LBFGS algorithm
        my_func = @(c_in) sharpness_figure_of_merit(c_in,Psi_in,Z,FOM);
        opt.max_objective = @(c_in) my_func(c_in);
        % Convergence criteria of the optimizaton step
        opt.ftol_rel = 1e-2;        
        opt.xtol_rel = 1e-2;   
        opt.maxeval = 500; 
        opt.verbose = 0;
        opt.lower_bounds = -2*ones(n_order,1);
        opt.upper_bounds = 2*ones(n_order,1);
        % Initial guess
        c_init = c_in;
        % Run optimization
        [c_in,~] = nlopt_optimize(opt,c_init);
        
        % Update
        r_out_update = exp(1i*Z*c_in).*r_out;
        
        % Build Psi_out matrix
        fprintf("Building Psi_out. \n")
        Psi_out = single(zeros(Nx^2,Nx^2));
        parfor ii = 1:Nx^2
            k_out_x = kx(ii); k_out_y = ky(ii);
            % Convert inputs to spatial basis
            r_update_ii = reshape(r_out_update(:,ii),Nx,Nx);
            r_k_r_ii = ifft(ifftshift(ifft(ifftshift(r_update_ii,1),Nx,1),2),Nx,2);
            % Convert the ouputs to spatial basis
            psi = reshape(r_k_r_ii,Nx,Nx).*exp(-1i*k_out_x*X-1i*k_out_y*Y);
        
            Psi_out(:,ii) = psi(:);
        end
        
        % Optimize outputs
        fprintf("Optimizing outputs.\n")
        my_func = @(c_out) sharpness_figure_of_merit(c_out,Psi_out,Z,FOM);
        opt.max_objective = @(c_out) my_func(c_out);
        % Initial guess
        c_init = c_out;
        % Run optimization
        [c_out,~] = nlopt_optimize(opt,c_init);
        
        % Update
        r_in_update = exp(1i*Z*c_out).*r_in;
        
        % Check the FOM
        I = abs(Psi_out*exp(1i*Z*c_out)).^2;
        M = fom(I,FOM);
        
        delta_M = abs((M-M_prev)/M_prev);

    end
    
    psi_opt = Psi_out*exp(1i*Z*c_out);
    psi_opt = flipud(fliplr(reshape(psi_opt,Nx,Nx)));
    
    figure(3)
    subplot 122
    imagesc(fliplr(interp2(abs(psi_opt).^2,3,'spline')));
    axis image
    colormap('hot')
    set(gca,'Visible','off')
    
    phi_in = Z*c_in; phi_out = Z*c_out;
   
    figure(4)
    subplot 121
    imagesc(wrapToPi(reshape(phi_in,Nx,Nx)));
    axis image
    colormap(colorcet('C1'))
    caxis([-pi pi])
    subplot 122
    imagesc(wrapToPi(reshape(phi_out,Nx,Nx)));
    axis image
    colormap(colorcet('C1'))
    caxis([-pi pi])

end
