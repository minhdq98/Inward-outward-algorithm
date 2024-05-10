function [stitch_mask] = build_stitching_mask(d_stitch, x_center, y_center, x, y)
    Nx = length(x);
    x_weight = zeros(1,Nx);
    x_weight(x >= x_center & x <= x_center+d_stitch) = 1-abs(x(x >= x_center & x <= x_center+d_stitch)-x_center)/d_stitch;
    x_weight(x <= x_center & x >= x_center-d_stitch) = 1-abs(x(x <= x_center & x >= x_center-d_stitch)-x_center)/d_stitch;
    y_weight = zeros(Nx,1);
    y_weight(y >= y_center & y <= y_center+d_stitch) = 1-abs(y(y >= y_center & y <= y_center+d_stitch)-y_center)/d_stitch;
    y_weight(y <= y_center & y >= y_center-d_stitch) = 1-abs(y(y <= y_center & y >= y_center-d_stitch)-y_center)/d_stitch;
    
    stitch_mask = y_weight.*x_weight;
end