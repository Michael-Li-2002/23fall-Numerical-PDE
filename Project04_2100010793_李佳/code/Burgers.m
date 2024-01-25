function pde = Burgers()
%BURGERS Burgers Equation
% f: conservative flux         Df: Dervivative of f    
% lbd/rbd: Left/Right Dirichlet boundary condition
% initu: initial value of u    exactu: exact (weak) solution
pde = struct('f',@f,'Df',@Df,'lbd',@lbd,'rbd',@rbd,'initu',@initu,'exactu',@weaku);

function w = f(u)
    w = u.^2/2;
end

function Deriv = Df(u)
    Deriv = u;
end

function ul = lbd(p)
    ul = 1 * ones(size(p,1),1);
end

function ur = rbd(p)
    ur = 0 * ones(size(p,1),1);
end

function u0 = initu(p)
    u0 = -min(sign(p),0);
end

function ureal = weaku(p)
    ureal = -min(sign(p(:,2)-p(:,1)/2),0);
end
end