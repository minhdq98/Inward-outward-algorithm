function [r_scale] = spatial_basis_truncation(X,Y,x_zone_scale,y_zone_scale,r_zone,if_do_svd,svd_thres)
    
    Nx = sqrt(height(r_zone));
    
    x_zone_scale_min = min(x_zone_scale,[],'all'); y_zone_scale_min = min(y_zone_scale,[],'all');
    x_zone_scale_max = max(x_zone_scale,[],'all'); y_zone_scale_max = max(y_zone_scale,[],'all');
    
    % Convert inputs to spatial basis
    r_k_r = single(zeros(Nx^2,Nx^2));
    for ii = 1:Nx^2
        r_k_ii = flipud(fliplr(reshape(r_zone(ii,:),Nx,Nx)));
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
    % Filter the out-of-zone elements
    mask_output = zeros(Nx^2,1);
    mask_output(X >= x_zone_scale_min & X <= x_zone_scale_max & Y >= y_zone_scale_min & Y <= y_zone_scale_max) = 1;
    mask_input = mask_output.';
    r_r = mask_output.*r_r.*mask_input;
    
    if if_do_svd
        [r_r] = do_svd(r_r,svd_thres,"spatial");
    end
    
    % Convert the inputs back to angular basis
    r_r_k = single(zeros(Nx^2,Nx^2));
    for ii = 1:Nx^2
        r_r_ii = reshape(r_r(ii,:),Nx,Nx);
        % Convert to angular basis
        r_r_k_ii = flipud(fliplr(fftshift(fft(fftshift(fft(r_r_ii,Nx,1),1),Nx,2),2)));
        r_r_k(ii,:) = transpose(r_r_k_ii(:));
    end
    % Convert the outputs back to angular basis
    r_scale = single(zeros(Nx^2,Nx^2));
    for ii = 1:Nx^2
        r_r_k_ii = reshape(r_r_k(:,ii),Nx,Nx);
        % Convert to angular basis
        r_k_ii = fftshift(fft(fftshift(fft(r_r_k_ii,Nx,1),1),Nx,2),2);
        r_scale(:,ii) = r_k_ii(:);
    end
end