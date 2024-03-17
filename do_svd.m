function [r_svd] = do_svd(r,svd_thres,type)
    Nx = sqrt(height(r));
    if type == "angular"
        % Convert to real space
        % Convert inputs to spatial basis
        r_k_r = single(zeros(Nx^2,Nx^2));
        for ii = 1:Nx^2
            r_k_ii = flipud(fliplr(reshape(r(ii,:),Nx,Nx)));
            r_k_r_ii = ifft(ifftshift(ifft(ifftshift(r_k_ii,1),Nx,1),2),Nx,2);
            r_k_r(ii,:) = transpose(r_k_r_ii(:));
        end
        % Convert the outputs to spatial basis
        r_r = single(zeros(Nx^2,Nx^2)); % The matrix whose both input and output are in spatial basis
        for ii = 1:Nx^2
            r_k_r_ii = reshape(r_k_r(:,ii),Nx,Nx);
            % Convert to spatial basis
            r_r_ii = ifft(ifftshift(ifft(ifftshift(r_k_r_ii,1),Nx,1),2),Nx,2);
            r_r(:,ii) = r_r_ii(:);
        end

        % Do SVD
        [u,s,v] = svd(r_r);
        svd_list = diag(s); svd_max = max(s,[],'all');
        N_svd = length(svd_list(svd_list >= svd_thres*svd_max));
        r_r_svd = u(:,1:N_svd)*s(1:N_svd,1:N_svd)*v(:,1:N_svd)';

        % Convert back to k space
        % Convert the inputs back to angular basis
        r_r_k = single(zeros(Nx^2,Nx^2));
        for ii = 1:Nx^2
            r_r_ii = reshape(r_r_svd(ii,:),Nx,Nx);
            % Convert to angular basis
            r_r_k_ii = flipud(fliplr(fftshift(fft(fftshift(fft(r_r_ii,Nx,1),1),Nx,2),2)));
            r_r_k(ii,:) = transpose(r_r_k_ii(:));
        end
        % Convert the outputs back to angular basis
        r_svd = single(zeros(Nx^2,Nx^2));
        for ii = 1:Nx^2
            r_r_k_ii = reshape(r_r_k(:,ii),Nx,Nx);
            % Convert to angular basis
            r_k_ii = fftshift(fft(fftshift(fft(r_r_k_ii,Nx,1),1),Nx,2),2);
            r_svd(:,ii) = r_k_ii(:);
        end
    elseif type == "spatial"
        % Do SVD
        [u,s,v] = svd(r);
        svd_list = diag(s); svd_max = max(s,[],'all');
        N_svd = length(svd_list(svd_list >= svd_thres*svd_max));
        r_svd = u(:,1:N_svd)*s(1:N_svd,1:N_svd)*v(:,1:N_svd)';
    end
    figure(10)
    plot(svd_list(svd_list > 0.01*svd_max));
end