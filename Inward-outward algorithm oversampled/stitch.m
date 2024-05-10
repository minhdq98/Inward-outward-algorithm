function [I_total,x_shift,y_shift] = stitch(I_total,I_zone,x_center,y_center,x,y,d_stitch)
    
    % Implement the stitching mask
    [stitch_mask] = build_stitching_mask(d_stitch, x_center, y_center, x, y);
    I_zone = I_zone.*stitch_mask;

    I_total = I_total+I_zone;
        
end