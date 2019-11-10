function [B,J1] = BMatrixXIGA2Tips(Xi,Eta,e,enrich_node,tip_node,...
                             xCr,crTip,alpha,N,dRdxi,dRdeta)
%                         
% same functionality as xigaBmatrix but for cracks with two tips
%

global controlPts element

sctr = element(e,:);
nn   = length(sctr);

% compute the jacobian of physical and parameter domain mapping
% then the derivative w.r.t spatial physical coordinates

pts        = controlPts(sctr,:);
jacob      = pts'*[dRdxi' dRdeta'];
J1         = det(jacob);
invJacob   = inv(jacob);
dRdx       = [dRdxi' dRdeta'] * invJacob;

Gpt = N * pts;                  % GP in global coord, used

% Bfem is always computed

Bfem = zeros(3,2*nn);
Bfem(1,1:2:2*nn)  = dRdx(:,1)';
Bfem(2,2:2:2*nn)  = dRdx(:,2)';
Bfem(3,1:2:2*nn)  = dRdx(:,2)';
Bfem(3,2:2:2*nn)  = dRdx(:,1)';

% Switch between non-enriched and enriched elements
if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    B = Bfem;
else                               % Enriched elements
    Bxfem = [] ;
    % loop on nodes, check node is enriched ...
    for in = 1 : nn
        if ( enrich_node(sctr(in)) == 1)     % H(x) enriched node
            % Enrichment function, H(x) at global Gauss point
            dist = signed_distance(xCr,Gpt);
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            dist = signed_distance(xCr,controlPts(sctr(in),:));
            Hi   = heaviside(dist);
            
            % Bxfem at node "in"
            BI_enr = [dRdx(in,1)*(Hgp - Hi) 0 ;
                      0 dRdx(in,2)*(Hgp - Hi) ;
                      dRdx(in,2)*(Hgp - Hi) dRdx(in,1)*(Hgp - Hi)];
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
        elseif ( enrich_node(sctr(in)) == 2) % B(x) enriched node 
            if ( tip_node(sctr(in))==1 )
                xTip  = crTip(1,:);
                alp = alpha(1,1);
            else
                xTip  = crTip(2,:);
                alp = alpha(1,2);
            end
            
            QT    = [cos(alp) sin(alp); 
                    -sin(alp) cos(alp)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi | theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = branch(r,theta,alp);
            
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(sctr(in),:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi | theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta);
            
            % components of Benr matrix
            
            aa = dRdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1) ;
            bb = dRdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1) ;
            B1_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2) ;
            bb = dRdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2) ;
            B2_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3) ;
            bb = dRdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3) ;
            B3_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4) ;
            bb = dRdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4) ;
            B4_enr = [aa 0 ; 0 bb ; bb aa];
            
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

