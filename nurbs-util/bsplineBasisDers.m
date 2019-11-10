Nxi     = [];
Neta    = [];
dNdxi   = [];
dNdeta  = [];
dRdxi   = [];
dRdeta  = [];
N       = [];

% compute derivative of basis functions w.r.t parameter coord

for in=1:noFnsU
    [Ni,dNi]  = NURBSbasis (connU(in),p,Xi,uKnot,weights);
    Nxi       = [Nxi Ni];
    dNdxi     = [dNdxi dNi];
end

for in=1:noFnsV
    [Ni,dNi]  = NURBSbasis (connV(in),q,Eta,vKnot,weights);
    Neta      = [Neta Ni];
    dNdeta    = [dNdeta dNi];
end

% derivate of R=Nxi*Neta w.r.t xi and eta
% this is derivative of shape functions in FEM

for j=1:noFnsV
    for i=1:noFnsU
        dRdxi  = [dRdxi  dNdxi(i) * Neta(j)];
        dRdeta = [dRdeta Nxi(i)   * dNdeta(j)];
        N      = [N      Nxi(i)   * Neta(j)];
    end
end