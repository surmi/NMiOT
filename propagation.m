function uo = propagation(ui,z , lamb,n0, dx)
    [Ny,Nx] = size(ui);
    Ui  = fftshift(fft2(ui));
    
    dfx = 1/Nx/dx;
    dfy = 1/Ny/dx;
    fx = (-Nx/2:Nx/2-1)*dfx;
    fy = (-Ny/2:Ny/2-1)*dfy;
    [Fx,Fy] = meshgrid(fx,fy);
    fz = sqrt((n0/lamb).^2 - Fx.^2 - Fy.^2);
    H = exp(1i*2*pi*fz*z);
    
    Uo = Ui.*H;
    
    uo = ifft2(fftshift(Uo));
end