% This function builds the Zernike wavefront matrices
function [Z] =  build_zernike(kmax,k,rad_order)
    %% 1. Build the Zernike matrix for reconstruction
    
    % 1.1. Build the normalized polar k space
    n1 = length(k);
    % Convert k to normalized polar coordinate
    k_norm_pol = zeros(n1,2); % 1st column: radial coordinate, 2nd column: angular coordinate
    for ii = 1:n1
        % radial coordinate
        k_norm_pol(ii,1) = sqrt(k(ii,1).^2+k(ii,2).^2)/kmax;
        sine_theta = k(ii,2)/sqrt(k(ii,1).^2+k(ii,2).^2); % Sine of the angular coordinate
        % Now, find the angular coordinate
        if k(ii,1) >= 0 && (k(ii,1) ~= 0 || k(ii,2) ~= 0)
            theta = asin(sine_theta);
        elseif k(ii,1) < 0 && (k(ii,1) ~= 0 || k(ii,2) ~= 0)
            theta = pi-asin(sine_theta);
        elseif k(ii,1)==0 && k(ii,2) == 0
            theta = 0;
        end
        % after this step, theta is from -pi/2 to 3*pi/2. wrap it to pi
        theta = wrapToPi(theta);
        k_norm_pol(ii,2) = theta;
    end
    
    % 1.2. Build the Zernike matrix
    % 1.2.1. Assign the radial and angular orders for the Zernike
    % polynomials
    order = []; % Initialize a matrix that contain the radial and angular orders 
    for ii = 1:rad_order
        for jj = -ii:2:ii
            order = vertcat(order,[jj ii]);
        end
    end
    n_order = length(order); % The number of Zernike polynomials
    % 1.2.2. Build matrix
    Z = zeros(n1,n_order);
    for jj = 1:n1
        rho = k_norm_pol(jj,1);
        theta = k_norm_pol(jj,2);
        m = order(:,1);
        n = order(:,2);
        if rho <= 1
            Z(jj,:) = zernfun(n,m,rho,theta,'norm');
        elseif rho > 1
            Z(jj,:) = 0;
        end
    end
    Z(:,1:2) = 0; % Remove tip and tilt
end