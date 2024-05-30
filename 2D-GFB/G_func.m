function g = G_func(r, x, landa)
    
    g = -1j/4.0*besselh(0, 2.0*pi/landa*norm(r - x));

end

