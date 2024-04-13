clear
clc
close all

%% 0. Provide the input info
% Save path 
path = "/media/minh/My Passport/Chicken breast Apr 8/Data/";

% Reflection matrix
r = load(""+path+"r.mat").r;

% k space
k = load(""+path+"/k.mat").k;
kx = k(:,1); ky = k(:,2); 
k_max = 4.1525;%max(abs(k),[],'all');

r(kx.^2+ky.^2 > k_max^2,:) = 0;
r(:,kx.^2+ky.^2 > k_max^2) = 0;

% Real space
L = 51.2;
dx = 2*pi/(2*k_max); dk = 2*pi/L;
x = 0:dx:L; y = x; Nx = length(x); [X,Y] = meshgrid(x,y); X = single(X(:)); Y = single(Y(:)); 
x_min = min(x,[],'all'); y_min = x_min;  x_max = max(x,[],'all'); y_max = x_max;

% FOM
FOM = 1;

% For inward optimization, how many inward steps do we want
n_inward = 12;

% Then, how many radial orders do we want for this inward optimization
rad_order_start = 18;
rad_order_step = 1;

% What is the shrinking criteria for the inward optimization?
FOM_thres = 0.95;

% Do we remove the previous zones in the outward optimization?
removal = 0;

% SVD threshold
svd_thres = 0.7;
svd_thres2 = 0.1;

%% 1. Build Zernike polynomials
Z_list = cell(1,1);

for step_id = 1:n_inward
    fprintf("Building Zernike matrix for the "+step_id+" inward step.\n")
    rad_order = rad_order_start+rad_order_step*(step_id-1);
    [Z] = build_zernike(k_max,k,rad_order);
    Z_list{step_id,1} = Z;       
end

%% 2. Inward optimization
I = zeros(Nx,Nx);
% Do SVD to remove the MS background
fprintf("Do SVD to remove the multiple scattering background, to improve the inward optimization.\n")
[r_svd] = do_svd(r,svd_thres,"angular");

r_step = r_svd;
phi_in_init = zeros(Nx^2,1); phi_out_init = zeros(Nx^2,1);
x_step = x; y_step = y;

for step_id = 1:n_inward

    Z = Z_list{step_id,1};
    
    % 2.1. Optimize the image at this step
    fprintf("Beging "+step_id+"th inward optimization.\n")
    [psi_step,phi_in_step,phi_out_step] = wavefront_correction(x,y,x_step,y_step,r_step,k,Z,FOM);
    fprintf("Done "+step_id+"th inward optimization.\n")
    
    % Update the reflection matrix
    r_step = exp(1i*phi_out_step).*r_step.*exp(1i*phi_in_step.');
    phi_in_init = phi_in_init+phi_in_step; phi_out_init = phi_out_init+phi_out_step;
            
    % Do the next steps if it is not the final step
    if step_id ~= n_inward
        % If this is step 1, located the sharpest spot on the FOV. 
        I_opt = abs(psi_step(:)).^2; I_max = max(I_opt,[],'all');
        if step_id == 1
            [ii_sharp] = find(I_opt == I_max);
            x_sharp = X(ii_sharp); y_sharp = Y(ii_sharp);
        end
        M_prev = fom(I_opt,FOM);

        % 2.2. Determine the size of the area surrounding the sharpest spot,
        % where FOM is lower than FOM_thres*M_prev
        fprintf("Beging finding IP size for the "+step_id+"th inward optimization.\n")
        ii = 1;
        x_step_min = max(0,x_sharp-ii*dx); y_step_min = max(0,y_sharp-ii*dx);
        x_step_max = min(L,x_sharp+ii*dx); y_step_max = min(L,y_sharp+ii*dx);
        I_step = I_opt(X >= x_step_min & X <= x_step_max & Y >= y_step_min & Y <= y_step_max);
        M_step = fom(I_step,FOM);      
        while M_step < FOM_thres*M_prev

            ii = ii+1;

            x_step_min = max(x_min,x_sharp-ii*dx); y_step_min = max(y_min,y_sharp-ii*dx);
            x_step_max = min(x_max,x_sharp+ii*dx); y_step_max = min(y_max,y_sharp+ii*dx);

            I_step = I_opt(X >= x_step_min & X <= x_step_max & Y >= y_step_min & Y <= y_step_max);

            M_step = fom(I_step,FOM);
        end
        fprintf("Done finding IP size for the "+step_id+"th inward optimization.\n")

        % 2.3. Now with the area determined, do spatial basis truncation
        x_step = x_step_min:dx:x_step_max; y_step = y_step_min:dx:y_step_max;
        fprintf("Beging spatial basis truncation for the "+step_id+"th inward optimization.\n")
        [r_step] = spatial_basis_truncation(X,Y,x_step,y_step,r_step,0,svd_thres);
        fprintf("Done spatial basis truncation for the "+step_id+"th inward optimization.\n")
    end
end

%% See how the FOV would look like if we apply the inward phase to the entire FOV
[r_svd2] = do_svd(r,svd_thres2,"angular");

% If for some reasons, the optimization shift the image transversely, manually reshift the image back
x_shift = 0; y_shift = 0; 
n = 0; % Compared to the final image size in the inward step, the progression step is how many pixels smaller?
d = 4*dx; % Compared to the progression step, the outward optimization zone is how much bigger?
x_step_min_shift = x_step_min-x_shift+n*dx; x_step_max_shift = x_step_max-x_shift-n*dx;
y_step_min_shift = y_step_min-y_shift+n*dx; y_step_max_shift = y_step_max-y_shift-n*dx;
r_shift = r_svd2.*exp(1i*(-kx.')*x_shift+1i*(-ky.')*y_shift);
r_shift = exp(1i*kx*x_shift+1i*ky*y_shift).*r_shift;

r_test = exp(1i*phi_out_init).*r_shift.*exp(1i*phi_in_init.');

psi_test = single(zeros(Nx,Nx));
% Find non-zero columns
[~,nz_in_ind] = find(r_test ~= 0); nz_in_ind = unique(nz_in_ind);
for iii = 1:length(nz_in_ind)
    ii = nz_in_ind(iii);
    k_in_x = kx(ii); k_in_y = ky(ii);
    % Convert outputs to spatial basis
    r_test_ii = reshape(r_test(:,ii),Nx,Nx);
    r_r_k_ii = ifft(ifftshift(ifft(ifftshift(r_test_ii,1),Nx,1),2),Nx,2);
    % Convert the inputs to spatial basis
    psi = r_r_k_ii.*reshape(exp(-1i*k_in_x*X-1i*k_in_y*Y),Nx,Nx);

    psi_test = psi_test+psi;
end
I_test = abs(psi_test).^2;
I_test = reshape(I_test,Nx,Nx);

% Build mask
x0 = 23; y0 = 21; % Central of the mask
wx0 = 15; wy0 = 15; % Width
rho = 0;
[gauss_mask] = build_mask(x,y,wx0,wy0,x0,y0,rho);
gauss_mask = reshape(gauss_mask,Nx,Nx);

figure(2)
imagesc(fliplr(interp2(I_test.*gauss_mask.^1,3,'spline')))
axis image
colormap('hot')

% Extract the image of the final inward optimization step
mask = zeros(Nx^2,1);
mask(X >= x_step_min_shift & X <= x_step_max_shift & Y >= y_step_min_shift & Y <= y_step_max_shift) = 1;
psi = psi_test.*reshape(mask,Nx,Nx);

% 3. Prepare for the outward optimization
% 3.1. Now the size of the progression step will be determined as
progression_step_by_2 = (x_step_max_shift-x_step_min_shift);

% The stitching mask parameter
d_stitch = d+progression_step_by_2/2;

% build the stitching mask for the inward zone
[stitch_mask] = build_stitching_mask(d_stitch, X(ii_sharp), Y(ii_sharp), x, y);

% 3.2. Number of zone in the first outward step
n_zone = 4;

% 3.3. The correction phase of the inward step
phi_prev = cell(1,2);
phi_prev{1,1} = phi_in_init; phi_prev{1,2} = phi_out_init;

% 3.4. The center of the zone in the next outward optimization, listed in
% clockwise order from min/min corner 
center_list = ones(n_zone,2);
center_list(1,:) = round([x_step_min_shift y_step_min_shift]/dx)*dx;
center_list(2,:) = round([x_step_min_shift y_step_max_shift]/dx)*dx;
center_list(3,:) = round([x_step_max_shift y_step_max_shift]/dx)*dx;
center_list(4,:) = round([x_step_max_shift y_step_min_shift]/dx)*dx;

% 3.5. What is the previous zone of these new 4 zones? All 1.
zone_prev_list = ones(4,2);

% 3.6. Store the coordinate of the previous zone
x_prev = cell(1,1); y_prev = cell(1,1);
x_prev{1,1} = x_step; y_prev{1,1} = y_step;

% 3.7. Zoning for the first outward optimization3
% Coordinates of each zones
list_x_zone = cell(1,1); list_y_zone = cell(1,1);
for zone_id = 1:n_zone
    x_center = center_list(zone_id,1); y_center = center_list(zone_id,2);
    
    x_zone = max(x_center-progression_step_by_2/2,x_min):dx:min(x_center+progression_step_by_2/2,x_max);
    y_zone = max(y_center-progression_step_by_2/2,y_min):dx:min(y_center+progression_step_by_2/2,y_max);
        
    list_x_zone{zone_id,1} = x_zone; list_y_zone{zone_id,1} = y_zone;
end

% Start counting the outward step from here
out_step = 0;

% Initialize the total image
I_total = I_test.*stitch_mask;

%% 4.Outward optimization, just run it until we are satisfied with the image
while out_step < 10 % Run until the image is good enough
    out_step = out_step+1;

    % 4.1. Optimizing each zone
    phi_zone_list = cell(1,1);
    
    FOM = 1;
    if_do_svd = 0; 
    svd_thres = 0.3;
    for zone_id = 1:n_zone
        % 4.1.1. Zone coordinate before scaling up
        x_zone = list_x_zone{zone_id,1}; y_zone = list_y_zone{zone_id,1};
        x_zone_min = min(x_zone,[],'all'); x_zone_max = max(x_zone,[],'all'); y_zone_min = min(y_zone,[],'all'); y_zone_max = max(y_zone,[],'all');
        x_center = center_list(zone_id,1); y_center = center_list(zone_id,2);
        
        % 4.1.2. Initial guess
        zone_prev_1 = zone_prev_list(zone_id,1); zone_prev_2 = zone_prev_list(zone_id,2);
        phi_in_zone_init = (phi_prev{zone_prev_1,1}+phi_prev{zone_prev_2,1})/2; 
        phi_out_zone_init = (phi_prev{zone_prev_1,2}+phi_prev{zone_prev_2,2})/2;
        phi_in_zone = phi_in_zone_init; phi_out_zone = phi_out_zone_init;

        % Update the initial guess
        r_zone = exp(1i*phi_out_zone_init).*r_shift.*exp(1i*phi_in_zone_init.');

        % 4.1.3. Optimize each zone
        % 4.1.3.1. Search for the zone size after scaling up for optimization
        x_zone_min_scale = max(x_zone_min-d,x_min);
        x_zone_max_scale = min(x_zone_max+d,x_max);
        x_zone_scale = x_zone_min_scale:dx:x_zone_max_scale;
        y_zone_min_scale = max(y_zone_min-d,y_min);
        y_zone_max_scale = min(y_zone_max+d,y_max);
        y_zone_scale = y_zone_min_scale:dx:y_zone_max_scale;

        % 4.1.3.2. Spatial basis truncation
        fprintf("Do spatial basis truncation before optimizing zone "+zone_id+" out of "+n_zone+" zones in the "+out_step+" outward step.\n")
        [r_zone] = spatial_basis_truncation(X,Y,x_zone_scale,y_zone_scale,r_zone,if_do_svd,svd_thres);
        fprintf("Done spatial basis truncation before optimizing zone "+zone_id+" out of "+n_zone+" zones in the "+out_step+" outward step.\n")

        % 4.1.3.3. Optimize the zone
        fprintf("Optimizing zone "+zone_id+" out of "+n_zone+" zones in the "+out_step+" outward step.\n")
        [psi_zone,phi_in_zone_step,phi_out_zone_step] = wavefront_correction(x,y,x_zone_scale,y_zone_scale,r_zone,k,Z,FOM);
        fprintf("Done optimizing zone "+zone_id+" out of "+n_zone+" zones in the "+out_step+" outward step.\n")

        % 4.1.3.4. Update the phase of the zone
        phi_in_zone = phi_in_zone+phi_in_zone_step; 
        phi_out_zone = phi_out_zone+phi_out_zone_step;
        r_zone = exp(1i*phi_out_zone_step).*r_zone.*exp(1i*phi_in_zone_step.');

        % If the previous zones are removed, rebuild the image with the
        % intact zone
        phi_zone_list{zone_id,1} = phi_in_zone; phi_zone_list{zone_id,2} = phi_out_zone;

        % 4.1.4. Remove the overlap area with the previous zones then put the 
        % zone image to the image of the big FOV
        I_zone = abs(psi_zone).^2;
        % build stitching mask
        [stitch_mask] = build_stitching_mask(d_stitch, x_center, y_center, x, y);
        I_total = I_total+I_zone.*stitch_mask;% 

        % Build mask
        x0 = 27; y0 = 20; % Central of the mask
        wx0 = 22; wy0 = 22; % Width
        rho = 0;
        [gauss_mask] = build_mask(x,y,wx0,wy0,x0,y0,rho);
        gauss_mask = reshape(gauss_mask,Nx,Nx);

        I_show = fliplr(interp2(I_total.*gauss_mask,3,'spline'));

        figure(1)
        imagesc(I_show)
        axis image
        colormap('hot')
        set(gca,'Visible','off')
        caxis([0 max(I_show,[],'all')])

    end

    % Update the range of the optimized region
    x_opt_min = x_max; x_opt_max = x_min; y_opt_min = y_max; y_opt_max = y_min;
    for zone_id = 1:n_zone
        x_zone = list_x_zone{zone_id,1}; y_zone = list_y_zone{zone_id,1};
        x_zone_min = min(x_zone,[],'all'); x_zone_max = max(x_zone,[],'all');
        y_zone_min = min(y_zone,[],'all'); y_zone_max = max(y_zone,[],'all');

        if x_zone_min < x_opt_min
            x_opt_min = x_zone_min;
        end
        if x_zone_max > x_opt_max
            x_opt_max = x_zone_max;
        end
        if y_zone_min < y_opt_min
            y_opt_min = y_zone_min;
        end
        if y_zone_max > y_opt_max
            y_opt_max = y_zone_max;
        end
    end
    x_opt_max = x_opt_max-1; x_opt_min = x_opt_min+1; y_opt_max = y_opt_max-1; y_opt_min = y_opt_min+1;

    % 4.2. Find the center of the zones in the next outward optimization step
    % along with their previous zones

    % 4.2.1. Let's first scan the old zones, but only those that have their centers
    % inside the FOV
    % These two table will store the center of the new zone along with their
    % previous zones 
    center_list_and_prev_zones_new = [];
    fprintf("Finding the centers of the new zones in the next outward optimization step.\n")

    for zone_id = 1:n_zone

        % 4.2.1.1. Here is the center of the old zone
        x_center = center_list(zone_id,1); y_center = center_list(zone_id,2);
        x_zone = list_x_zone{zone_id,1}; y_zone = list_y_zone{zone_id,1};
        x_zone_min = min(x_zone,[],'all'); y_zone_min = min(y_zone,[],'all');
        x_zone_max = max(x_zone,[],'all'); y_zone_max = max(y_zone,[],'all');
        lx_zone = x_zone_max-x_zone_min; ly_zone = y_zone_max-y_zone_min;

        % Only do the next step if the center of this old zone is inside the
        % FOV
        if x_center >= x_min & x_center <= x_max & y_center >= y_min & y_center <= y_max

            % 4.2.1.2. The center of the new zones spawning from. These new center points
            % will belong to 2 or 3 of the following 4 points. Note that the center
            % point here can be out of the FOV.
            for ii = 1:4
                switch ii
                    case 1
                        x_center_new = x_center-round(progression_step_by_2/2/dx)*dx;
                        y_center_new = y_center-round(progression_step_by_2/2/dx)*dx;
                    case 2
                        x_center_new = x_center+round(progression_step_by_2/2/dx)*dx;
                        y_center_new = y_center-round(progression_step_by_2/2/dx)*dx;
                    case 3
                        x_center_new = x_center+round(progression_step_by_2/2/dx)*dx;
                        y_center_new = y_center+round(progression_step_by_2/2/dx)*dx;
                    case 4
                        x_center_new = x_center-round(progression_step_by_2/2/dx)*dx;
                        y_center_new = y_center+round(progression_step_by_2/2/dx)*dx;
                end

                % 4.2.1.3. Then, check if the the new center is in the outward direction. If
                % the new center belongs to the older zone then it's not the
                % outward direction
                if (x_center_new < x_opt_min | x_center_new > x_opt_max) | (y_center_new < y_opt_min | y_center_new > y_opt_max)
                    x_zone_min_new = max(min(x_center_new:-dx:x_center_new-progression_step_by_2/2),x_min);
                    x_zone_max_new = min(max(x_center_new:dx:x_center_new+progression_step_by_2/2),x_max);
                    y_zone_min_new = max(min(y_center_new:-dx:y_center_new-progression_step_by_2/2),y_min);
                    y_zone_max_new = min(max(y_center_new:dx:y_center_new+progression_step_by_2/2),y_max);

                    lx_zone_new = x_zone_max_new-x_zone_min_new; ly_zone_new = y_zone_max_new-y_zone_min_new;

                    if x_center_new > 0 & x_center_new < L & y_center_new > 0 & y_center_new < L
                        center_list_and_prev_zones_new = vertcat(center_list_and_prev_zones_new,...
                                        [x_center_new y_center_new zone_id]);
                    end
                end
            end
        end

    end

    % 4.2.2. After this step, we get a list of center points of the new zones and
    % their previous zones. Here is the list of center points only
    center_list_new = center_list_and_prev_zones_new(:,1:2);

    % 4.2.3. The new list of centers becomes the current list of centers
    center_list = uniquetol(round(center_list_new/dx),0.05,'ByRows',true)*dx;

    % 4.2.4. The number of new zones is temporarily
    n_zone = length(center_list);

    % 4.2.5. The current zones now become the previous zones
    x_prev = list_x_zone; y_prev = list_y_zone;

    % 4.2.6. The current list of correction phases becomes the previous list of
    % current phase
    phi_prev = phi_zone_list;

    % 4.2.7. Now scan through the new zones, find their respective coordinates and
    % their two previous zones
    fprintf("Finding the coordinates and previous zones of the new zones in the next outward optimization step.\n")
    list_x_zone = cell(1,1); list_y_zone = cell(1,1); zone_prev_list = zeros(n_zone,2);
    for zone_id = 1:n_zone
        x_center = center_list(zone_id,1); y_center = center_list(zone_id,2);
        x_zone_min = max(min(x_center:-dx:x_center-progression_step_by_2/2),x_min);
        x_zone_max = min(max(x_center:dx:x_center+progression_step_by_2/2),x_max);
        y_zone_min = max(min(y_center:-dx:y_center-progression_step_by_2/2),y_min);
        y_zone_max = min(max(y_center:dx:y_center+progression_step_by_2/2),y_max);
        x_zone = x_zone_min:dx:x_zone_max; y_zone = y_zone_min:dx:y_zone_max;
        list_x_zone{zone_id,1} = x_zone; list_y_zone{zone_id,1} = y_zone; 
    end
    % Now find the previous zones
    for zone_id = 1:n_zone
        jj = 0;
        for ii = 1:length(center_list_and_prev_zones_new)

            if round(center_list_and_prev_zones_new(ii,1)/dx) == round(center_list(zone_id,1)/dx) & ...
                    round(center_list_and_prev_zones_new(ii,2)/dx) == round(center_list(zone_id,2)/dx)
                jj = jj+1;
                zone_prev_list(zone_id,jj) = center_list_and_prev_zones_new(ii,3);
            end
            % If this zone only takes one previous zone
            if jj == 1
                zone_prev_list(zone_id,2) = zone_prev_list(zone_id,1);
            end
        end
    end

    save(""+path+"I_io.mat",'I_total')
    
    if n_zone == 0
        break
    end
end
