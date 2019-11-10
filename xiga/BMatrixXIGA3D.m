function [B,J1] = BMatrixXIGA3D(e,enrich_node,...
                                N, dNdxi, dNdeta, dNdzeta,pt)
% Compute the strain-displacement B matrix for 3D enriched IGA
% B = [Bfem | Benr]
% Input:
%   e: element under consideration
%   enrich_node: enriched nodes
%   N: NURBS shape functions 
%   dNdxi: derivatives of N w.r.t xi
%   pt: (xi,eta,zeta) in the B8 parent element
% 
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

global controlPts element elementV xCr xTip levelSets levelSetsB8 node

sctr = element(e,:);
nn   = length(sctr);

% compute the jacobian of physical and parameter domain mapping
% then the derivative w.r.t spatial physical coordinates

pts        = controlPts(sctr,:);
jacob      = pts'*[dNdxi' dNdeta' dNdzeta'];
J1         = det(jacob);
invJacob   = inv(jacob);
dNdx       = [dNdxi' dNdeta' dNdzeta'] * invJacob;
Gpt        = N * pts;                  % GP in global coord, used
Phi        = levelSets(sctr,:);

% B8 element for level sets

% [N8,dN8dxi] = lagrange_basis('B8',pt); 
% sctrB8      = elementV(e,:);
% J0          = node(sctrB8,:)'*dN8dxi;
% invJ0       = inv(J0);
% dN8dx       = dN8dxi*invJ0;
% Phi         = levelSetsB8(sctrB8,:);

% Bfem is always computed

nn3 = 3*nn;

Bfem = zeros(6,nn3);
Bfem(1,1:3:nn3)  = dNdx(:,1)';
Bfem(2,2:3:nn3)  = dNdx(:,2)';
Bfem(3,3:3:nn3)  = dNdx(:,3)';

Bfem(4,1:3:nn3)  = dNdx(:,2)';
Bfem(4,2:3:nn3)  = dNdx(:,1)';

Bfem(5,2:3:nn3)  = dNdx(:,3)';
Bfem(5,3:3:nn3)  = dNdx(:,2)';

Bfem(6,1:3:nn3)  = dNdx(:,3)';
Bfem(6,3:3:nn3)  = dNdx(:,1)';

% Switch between non-enriched and enriched elements

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    B = Bfem;
else                               % Enriched elements
    Bxfem = [] ;
    % loop on nodes, check node is enriched ...
    
    for in = 1 : nn
        nodeId  = sctr(in);
        enrNoId = enrich_node(nodeId);
        
        if ( enrNoId == 1)     % H(x) enriched node
            % Enrichment function, H(x) at global Gauss point
            dist = Gpt(1,3)-xCr(1,2);
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            dist = controlPts(nodeId,3)-xCr(1,2);
            Hi   = heaviside(dist); Hi=0;
            HH   = Hgp - Hi;
            aa   = dNdx(in,1)*HH;
            bb   = dNdx(in,2)*HH;
            cc   = dNdx(in,3)*HH;
            
            % Bxfem at node "in"
            
            BI_enr = [aa 0  0;
                      0  bb 0;
                      0  0  cc;
                      bb aa 0;
                      0 cc bb;
                      cc 0 aa];
            
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
        elseif ( enrNoId == 2) % B(x) enriched node   
            % branch functions in terms of B8 level sets
            
%             [Br,dBdxi,dBdeta,dBdzeta] = branch3d(Phi,N8,...
%                         dN8dxi(:,1),dN8dxi(:,2),dN8dxi(:,3));
%             
%             dBdx  = [dBdxi' dBdeta' dBdzeta'] * invJ0;
            
            % branch functions in terms of NURBS level sets

            [Br,dBdxi,dBdeta,dBdzeta] = branch3d(Phi,N',...
                                      dNdxi',dNdeta',dNdzeta');
            
            dBdx  = [dBdxi' dBdeta' dBdzeta'] * invJacob;
            
            % compute branch functions at node "in"
            
            %[BrI] = branchNode(levelSets(nodeId,1),levelSets(nodeId,2));
            BrI = [0 0 0 0];
            % components of Benr matrix
            
            Nin = N(in);
            
            % first branch function
            
            aa = dNdx(in,1)*(Br(1)-BrI(1)) + Nin*dBdx(1,1);
            bb = dNdx(in,2)*(Br(1)-BrI(1)) + Nin*dBdx(1,2);
            cc = dNdx(in,3)*(Br(1)-BrI(1)) + Nin*dBdx(1,3);
            
            B1_enr = [aa 0 0;
                      0 bb 0;
                      0 0 cc;
                      bb aa 0;
                      0 cc bb;
                      cc 0 aa];
            
            % second branch function
            
            aa = dNdx(in,1)*(Br(2)-BrI(2)) + Nin*dBdx(2,1);
            bb = dNdx(in,2)*(Br(2)-BrI(2)) + Nin*dBdx(2,2);
            cc = dNdx(in,3)*(Br(2)-BrI(2)) + Nin*dBdx(2,3);
            
            B2_enr = [aa 0 0;
                      0 bb 0;
                      0 0 cc;
                      bb aa 0;
                      0 cc bb;
                      cc 0 aa];
            
            % third branch function
            
            aa = dNdx(in,1)*(Br(3)-BrI(3)) + Nin*dBdx(3,1);
            bb = dNdx(in,2)*(Br(3)-BrI(3)) + Nin*dBdx(3,2);
            cc = dNdx(in,3)*(Br(3)-BrI(3)) + Nin*dBdx(3,3);
            
            B3_enr = [aa 0 0;
                      0 bb 0;
                      0 0 cc;
                      bb aa 0;
                      0 cc bb;
                      cc 0 aa];
            
            % fourth branch function
            
            aa = dNdx(in,1)*(Br(4)-BrI(4)) + Nin*dBdx(4,1);
            bb = dNdx(in,2)*(Br(4)-BrI(4)) + Nin*dBdx(4,2);
            cc = dNdx(in,3)*(Br(4)-BrI(4)) + Nin*dBdx(4,3);
            
            B4_enr = [aa 0 0;
                      0 bb 0;
                      0 0 cc;
                      bb aa 0;
                      0 cc bb;
                      cc 0 aa];
            
            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
        end
    end          % end of loop on nodes
    % B matrix
    B = [Bfem Bxfem];
    clear Bfem; clear Bxfem;
end              % end of switch between enriched and non-enriched elements

