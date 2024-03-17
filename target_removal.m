function [r] = target_removal(r,mask_remove)
    
    Nx = sqrt(height(r));
    
    % Find non-zero rows
    [nz_row_ind,~] = find(r ~= 0); nz_row_ind = unique(nz_row_ind); 
    % Convert inputs to spatial basis
    r_k_r = single(zeros(Nx^2,Nx^2));
    for iii = 1:length(nz_row_ind)
        ii = nz_row_ind(iii);
        r_k_ii = flipud(fliplr(reshape(r(ii,:),Nx,Nx)));
        r_k_r_ii = ifft(ifftshift(ifft(ifftshift(r_k_ii,1),Nx,1),2),Nx,2);
        r_k_r(ii,:) = transpose(r_k_r_ii(:));
    end
    % Find non-zero columns
    [~,nz_col_ind] = find(r_k_r ~= 0); nz_col_ind = unique(nz_col_ind);
    % Convert the outputs to spatial basis
    r_r = single(zeros(Nx^2,Nx^2)); % The matrix whose both input and output are in spatial basis
    for iii = 1:length(nz_col_ind)
        ii = nz_col_ind(iii);
        r_k_r_ii = reshape(r_k_r(:,ii),Nx,Nx);
        % Convert to spatial basis
        r_r_ii = ifft(ifftshift(ifft(ifftshift(r_k_r_ii,1),Nx,1),2),Nx,2);
        r_r(:,ii) = r_r_ii(:);
    end

    % Remove
    r_r = r_r.*mask_remove;
    
    % Find non-zero rows
    [nz_row_ind,~] = find(r_r ~= 0); nz_row_ind = unique(nz_row_ind); 
    % Convert the inputs back to angular basis
    r_r_k = single(zeros(Nx^2,Nx^2));
    for iii = 1:length(nz_row_ind)
        ii = nz_row_ind(iii);
        r_r_ii = reshape(r_r(ii,:),Nx,Nx);
        % Convert to angular basis
        r_r_k_ii = flipud(fliplr(fftshift(fft(fftshift(fft(r_r_ii,Nx,1),1),Nx,2),2)));
        r_r_k(ii,:) = transpose(r_r_k_ii(:));
    end
    % Find non-zero columns
    [~,nz_col_ind] = find(r_r_k ~= 0); nz_col_ind = unique(nz_col_ind);
    % Convert the outputs back to angular basis
    r = single(zeros(Nx^2,Nx^2));
    for iii = 1:length(nz_col_ind)
        ii = nz_col_ind(iii);
        r_r_k_ii = reshape(r_r_k(:,ii),Nx,Nx);
        % Convert to angular basis
        r_k_ii = fftshift(fft(fftshift(fft(r_r_k_ii,Nx,1),1),Nx,2),2);
        r(:,ii) = r_k_ii(:);
    end
end