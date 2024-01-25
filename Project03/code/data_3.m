function pde = data_3
%DATA_2 : initial data is only piecewise continuous.
pde = struct('initdata',@initdata,'exactu',@exactu);

function u0 = initdata(p)
    u0 = 1/2* p .*(1 + sign(1/2-p));
end

function u = exactu(p)
    trun = 3;
    serie = (1:trun);
    u = sum((-cos((serie*pi)/2)./(serie*pi)+2*sin((serie*pi)/2)./(serie.^2*pi^2))...
        .*sin(pi*p(:,2)*serie).*exp(-pi^2*p(:,1)*serie.^2),2);
end
end



