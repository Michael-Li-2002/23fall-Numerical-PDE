function pde = data_2
%DATA_2 : initial data continuous, but its derivative is not.
pde = struct('initdata',@initdata,'exactu',@exactu);

function u0 = initdata(p)
    u0 = sin(2*pi*min(p,1/2));
end

function u = exactu(p)
    trun = 3;
    serie = (1:trun);
    u = 1/2 * exp(-4*pi^2*p(:,1)).*sin(2*pi*p(:,2))...
        + sum(4*(-1).^serie./(pi*(2*serie+1).*(2*serie-3))...
        .* sin(p(:,2)*(2*serie-1)*pi).* exp(-p(:,1)*(2*serie-1).^2*pi^2),2);
end
end

