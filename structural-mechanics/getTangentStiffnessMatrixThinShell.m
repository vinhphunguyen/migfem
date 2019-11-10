function [K,fint] = getTangentStiffnessMatrixThinShell (U)
% compute the tangent stiffness matrix and the internal force
% vector for Kirchhoff-Love Shell Total Lagrangian formulation for finite deformation.
% VP Nguyen
% Cardiff University, March 2013.
% Todo:
% 1. COnvert to Bezier extraction: DONE
% 2. Clean the code
% 3. Is it better to replace Newton-Raphson by explicit dynamics combined
% with dynamic relaxation?

global element p q index elRangeU elRangeV elRangeW uKnot vKnot wKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP nu E memStiff benStiff K

global shapes  gradsx  gradse  grads2x  grads2e  grads2xe

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*3;
nElmLK = nElDof^2;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

noElems=size(element,1);
nSprGK = nElmLK*noElems;
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);


fint = zeros(noDofs,1);        % internal force vector

noBasis = nElNod;

% Loop over elements (knot spans)

for e=1:noElems
    sctr   = element(e,:);          %  element connectivity
    nn     = length(sctr);
    
    nn3      = 3*nn;
    sctrB = zeros(1,nn3);
    
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr;
    
    sctrB(1:3:nn3) = sctrx;
    sctrB(2:3:nn3) = sctry;
    sctrB(3:3:nn3) = sctrz;
        
    elemDisp = [U(sctrx) U(sctry) U(sctrz)];
    
    pts    = controlPts(sctr,:);   % initial coords of elemenrt nodes
    ptsc   = pts + elemDisp;       % current coords of element nodes
           
    Ke = Ke0;
        
    for gp=1:size(W,1)      % loop over Gauss points
        wt      = W(gp);
        
        R       = shapes  {e}(gp,:);
        dRdxi   = gradsx  {e}(gp,:);
        dRdeta  = gradse  {e}(gp,:);
        dR2dxi  = grads2x {e}(gp,:);
        dR2det  = grads2e {e}(gp,:);
        dR2dxe  = grads2xe{e}(gp,:);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob   = [dRdxi; dRdeta]          * pts; % 2x3 matrix
        jacob2  = [dR2dxi; dR2det; dR2dxe] * pts; % 3x3 matrix
        jacobc  = [dRdxi; dRdeta]          * ptsc; % 2x3 matrix
        jacob2c = [dR2dxi; dR2det; dR2dxe] * ptsc; % 3x3 matrix
        
        % a1, a2 and a3 vectors (surface basis vectors)
        % and its derivatives
        
        a1b    = jacob(1,:);
        a2b    = jacob(2,:);
        a3b    = cross(a1b,a2b);
        normb  = norm(a3b);
        a3b    = a3b/normb; J1    = normb;
        
        a11b   = jacob2(1,:);
        a22b   = jacob2(2,:);
        a12b   = jacob2(3,:);
        
        a11   = jacob2c(1,:);
        a22   = jacob2c(2,:);
        a12   = jacob2c(3,:);
        
        % basis vectors in deformed configuration
        
        a1 = jacobc(1,:);
        a2 = jacobc(2,:);
        a3 = cross(a1,a2);
        norma = norm(a3);
        a3    = a3/norma;
        normi = 1/norma;
        
        % compute the constitutive matrix C
        % reference configuration
        
        a_11 = dot(a1b,a1b); a_12 = dot(a1b,a2b);
        a_21 = dot(a2b,a1b); a_22 = dot(a2b,a2b);
        
        aa1 = [a_11 a_21; a_12 a_22 ] \ [1;0];
        aa2 = [a_11 a_21; a_12 a_22 ] \ [0;1];
        
        au11 = aa1(1);
        au12 = aa1(2);
        au22 = aa2(2);
        
        C = [au11^2 nu*au11*au22+(1-nu)*au12^2 au11*au12;
            nu*au11*au22+(1-nu)*au12^2 au22^2 au22*au12;
            au11*au12 au22*au12 0.5*((1-nu)*au11*au22+(1+nu)*au12^2)];
        
        % compute the membrane and bending stress resultant
        % first, membrane and bending strains
        
        alpha=0.5*[   dot(a1,a1) - dot(a1b,a1b);
            dot(a2,a2) - dot(a2b,a2b);
            2*(dot(a1,a2) - dot(a1b,a2b))];
        
        beta = [   dot(a11b,a3b) - dot(a11,a3);
            dot(a22b,a3b) - dot(a22,a3);
            2*(dot(a12b,a3b) - dot(a12,a3))];
        
        n = memStiff * C * alpha;
        m = benStiff * C * beta;
        
        % membrane and bending B matrices
        
        da3a11  = dot(a3,a11);
        da3a22  = dot(a3,a22);
        da3a12  = dot(a3,a12);
        
        a11ca2  = cross(a11,a2);
        a1ca11  = cross(a1, a11);
        a22ca2  = cross(a22,a2);
        a1ca22  = cross(a1,a22);
        a3ca1   = cross(a3,a1);
        a2ca3   = cross(a2,a3);
        a12ca2  = cross(a12,a2);
        a1ca12  = cross(a1,a12);
        
        % weitghts
        
        wi = J1 * wt;
        
        Hphi   = [a11;a22;2*a12];
        mthphi = m'*Hphi;
        
        tx   = [0 -a3(3) a3(2);
            a3(3) 0 -a3(1);
            -a3(2) a3(1) 0];
        
        mthphix   = [0 -mthphi(3) mthphi(2);
            mthphi(3) 0 -mthphi(1);
            -mthphi(2) mthphi(1) 0];
                                
        for i =1:noBasis
            sctri =3*i-2 : 3*i;
            dRIdx = dRdxi (i);
            dRIdy = dRdeta(i);
            
            wa    = -dRIdx*a2 + dRIdy*a1;            
            wax   = [0 -wa(3) wa(2);
                wa(3) 0 -wa(1);
                -wa(2) wa(1) 0];
            
            %Hpa   = [dR2dxi(i) dR2det(i) 2*dR2dxe(i)];
            
            dotmHpa=m(1)*dR2dxi(i) + m(2)*dR2det(i) + m(3)*2*dR2dxe(i);
            
            dJdpI = dRIdx * a2ca3 + dRIdy * a3ca1;
                        
            BmemI=[dRIdx*a1;
                dRIdy*a2;
                dRIdy*a1 + dRIdx*a2];
            
            BI1 = -dR2dxi(i)*a3 + normi*(dRIdx*a11ca2 + dRIdy*a1ca11 + ...
                da3a11*(dRIdx*a2ca3 + dRIdy*a3ca1));
            
            BI2 = -dR2det(i)*a3 + normi*(dRIdx*a22ca2 + dRIdy*a1ca22 + ...
                da3a22*(dRIdx*a2ca3 + dRIdy*a3ca1));
            
            BI3 = -dR2dxe(i)*a3 + normi*(dRIdx*a12ca2 + dRIdy*a1ca12 + ...
                da3a12*(dRIdx*a2ca3 + dRIdy*a3ca1));
            
            BbenI=[BI1;BI2;2*BI3];
            
            sctrI = [3*sctr(i)-2 3*sctr(i)-1 3*sctr(i)] ;
            fint(sctrI) = fint(sctrI)  +  BmemI' * n * wi + BbenI' * m * wi;
            
            for j=1:noBasis
                dRJdx = dRdxi (j);
                dRJdy = dRdeta(j);
                sctrj =3*j-2 : 3*j;
                
                %Hpb   = [dR2dxi(j) dR2det(j) 2*dR2dxe(j)];
                
                dotmHpb=m(1)*dR2dxi(j) + m(2)*dR2det(j) + m(3)*2*dR2dxe(j);
                
                aa = dRIdx*dRJdx*n(1) + dRIdy*dRJdy*n(2) + ...
                    (dRIdx*dRJdy + dRJdx*dRIdy )*n(3);
                
                %aa = dRIdRJ1(i,j) + dRIdRJ2(i,j) + dRIdRJ3(i,j);
                
                dJdpJ = dRJdx * a2ca3 + dRJdy * a3ca1;                
                wb    = -dRJdx*a2 + dRJdy*a1;                
                wbx   = [ 0 -wb(3) wb(2);
                    wb(3) 0 -wb(1);
                    -wb(2) wb(1) 0];
                
                waxwbx=[wa(3)*wb(3)+wa(2)*wb(2) -wa(2)*wb(1) -wa(3)*wb(1);
                    -wa(1)*wb(2) wa(3)*wb(3)+wa(1)*wb(1) -wa(3)*wb(2);
                    -wa(1)*wb(3) -wa(2)*wb(3) wa(2)*wb(2)+wa(1)*wb(1)];
                
                BmemJ=[dRJdx*a1;dRJdy*a2;dRJdy*a1 + dRJdx*a2];
                
                BI1 = -dR2dxi(j)*a3 + normi*(dRJdx*a11ca2 + dRJdy*a1ca11 + ...
                    da3a11*(dRJdx*a2ca3 + dRJdy*a3ca1));
                
                BI2 = -dR2det(j)*a3 + normi*(dRJdx*a22ca2 + dRJdy*a1ca22 + ...
                    da3a22*(dRJdx*a2ca3 + dRJdy*a3ca1));
                
                BI3 = -dR2dxe(j)*a3 + normi*(dRJdx*a12ca2 + dRJdy*a1ca12 + ...
                    da3a12*(dRJdx*a2ca3 + dRJdy*a3ca1));
                
                BbenJ=[BI1;BI2;2*BI3];
                dd =  -dRIdx*dRJdy + dRIdy*dRJdx;
                m'*Hphi*a3'
                tem = - dJdpI'*(BbenJ*m)' - BbenI'*m * dJdpJ...
                    + m'*Hphi*a3'*(-normi*dJdpI*dJdpJ' + normi* waxwbx + dd*tx )...
                    - dotmHpa*wbx - dotmHpb*wax'-dd*mthphix;
                
                Ke(sctri,sctrj) =  Ke(sctri,sctrj) + aa*eye(3)*wi + tem*normi*wi ...
                                   + memStiff * BmemI' * C * BmemJ * wi + ...
                                     benStiff * BbenI' * C * BbenJ * wi;
            end
        end                                       
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position
    
end

% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);
