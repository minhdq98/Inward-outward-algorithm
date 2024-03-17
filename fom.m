function M = fom(I,FOM)
    if FOM == 'entropy'
        M = sum(I.*log(I),'all');
    else
        M = sum(I.^FOM,'all');
    end
end