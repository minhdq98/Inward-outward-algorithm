path = "/media/minh/WD_BLACK/Chicken breast march 15/Data/";

r = load("/media/minh/WD_BLACK/Chicken breast march 15/Data/r_z_2D_2.mat").r_z;
% k space
k = single(load("/media/minh/WD_BLACK/Chicken breast march 15/k.mat").k);
k_max = max(sqrt(k(:,1).^2+k(:,2).^2),[],'all'); N = length(k);
kx = k(:,1); ky = k(:,2);
N = length(k);

L = 51.2;
dx = 2*pi/(2*k_max);
x = dx/2:dx:L; y = x; Nx = length(x); [X,Y] = meshgrid(x,y); X = single(X(:)); Y = single(Y(:)); 

% Create a square k space
dk = 2*pi/L;
kx_new = -k_max+dk/2:dk:k_max; ky_new = kx_new;
[Kx_new,Ky_new] = meshgrid(kx_new,ky_new);
kx_new = single(Kx_new(:)); ky_new = single(Ky_new(:)); 
N_new = length(kx_new);

%% Convert inputs to spatial basis
fprintf("Convert inputs to spatial basis.\n")
r_k_r = single(zeros(N,Nx^2));
for ii = 1:N
    r_k_r(ii,:) = finufft2d3(ky,kx,r(ii,:).',-1,1e-2,Y,X).'/Nx;
end
% Convert the outputs to spatial basis
fprintf("Convert the outputs to spatial basis.\n")
r_r = single(zeros(Nx^2,Nx^2)); % The matrix whose both input and output are in spatial basis
for ii = 1:Nx^2
    r_r(:,ii) = finufft2d3(-ky,-kx,r_k_r(:,ii),-1,1e-2,Y,X)/Nx;
end

%% Remove off-diagonal
w = 10;
for ii = 1:Nx^2
    x_in = X(ii); y_in = Y(ii);
    window = zeros(Nx^2,1);
    window((X-x_in).^2+(Y-y_in).^2 <= w^2) = 1;
    r_r(:,ii) = r_r(:,ii).*window;

end

% Remove reflection of the objective lens at the edge of the image
x_min = x(1); x_max = x(end); y_min = y(1); y_max = y(end);
d_remove = 0*dx;
r_r(X <= d_remove | X >= L-d_remove | Y <= d_remove | Y >= L-d_remove,:) = 0;
r_r(:,X <= d_remove | X >= L-d_remove | Y <= d_remove | Y >= L-d_remove) = 0;

% SVD
[u,s,v] = svd(r_r);
N_sv = 400;
r_r = u(:,1:N_sv)*s(1:N_sv,1:N_sv)*v(:,1:N_sv)';

%% Convert the inputs back to angular basis
fprintf("Convert the inputs back to angular basis.\n")
r_r_k = single(zeros(Nx^2,N_new));
for ii = 1:Nx^2
    r_r_k(ii,:) = finufft2d3(Y,X,r_r(ii,:).',1,1e-2,ky_new,kx_new)/Nx;
end
% Convert the outputs back to angular basis
fprintf("Convert the outputs back to angular basis.\n")
r = single(zeros(N_new,N_new));
for ii = 1:N
    r_ii = finufft2d3(Y,X,r_r_k(:,ii),1,1e-2,-ky_new,-kx_new)/Nx;
    
    % Everything outside NA is remove
    if kx_new(ii)^2+ky_new(ii)^2 > k_max^2
        r_ii = 0*r_ii;
    end
    
    r_ii(kx_new.^2+ky_new.^2 > k_max^2) = 0;
    
    r(:,ii) = r_ii(:);
    
end

% Save the results
k = [kx_new ky_new];
save(""+path+"r.mat",'r')
save(""+path+"k.mat",'k','dx')

%% Convert to old angular basis
fprintf("Convert the inputs back to angular basis.\n")
r_r_k = single(zeros(Nx^2,N));
for ii = 1:Nx^2
    r_r_k(ii,:) = finufft2d3(Y,X,r_r(ii,:).',1,1e-2,kyr_,kx)/Nx;
end
% Convert the outputs back to angular basis
fprintf("Convert the outputs back to angular basis.\n")
r = single(zeros(N,N));
for ii = 1:N
    r(:,ii) = finufft2d3(Y,X,r_r_k(:,ii),1,1e-2,-ky,-kx)/Nx;
     
end

% Save the results
k = [kx ky];
save(""+path+"r_circ.mat",'r')
save(""+path+"k_circ.mat",'k','dx')
