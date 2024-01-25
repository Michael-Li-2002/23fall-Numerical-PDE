function pde = data_1
%DATA_1 : smooth initial data
pde = struct('initdata',@initdata,'exactu',@exactu);

function u0 = initdata(p)
    u0 = sin(pi*p);
end

function u = exactu(p)
    u = exp(-pi^2*p(:,1)).*sin(pi*p(:,2));
end
end

