function [M,grad] = sharpness_figure_of_merit(c,Psi,Z,FOM)

    % Calculate the FOM
    c = c(:);
    phi = Z*c;
    psi = Psi*exp(1i*phi);
    I = abs(psi).^2;
   
    M = fom(I,FOM);
    
    phi = double(phi);
    
% Calculate the gradient
    if nargout>1
        if FOM == 'entropy'
            grad = -double(2*imag(((log(I)+1).*psi)'*Psi.*(exp(1i*phi.')))*Z);
        else
            grad = -double(2*FOM*imag((I.^(FOM-1).*psi)'*Psi.*(exp(1i*phi.')))*Z);
        end
    end
end
