function [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta)
% Compute the strain-displacement B matrix for enriched IGA elements
% Vinh Phu Nguyen, nvinhphu@gmail.com
% Delft University of Technology, The Netherlands
% 2012

global controlPts element xCrack xTips crack_node ...
       inclusion_node VOIDs CHI ep weakEnrFunc

sctr = element(e,:);
nn   = length(sctr);

% compute the jacobian of physical and parameter domain mapping
% then the derivative w.r.t spatial physical coordinates

pts        = controlPts(sctr,:);
jacob      = pts'*[dRdxi' dRdeta'];
J1         = det(jacob);
invJacob   = inv(jacob);
dRdx       = [dRdxi' dRdeta'] * invJacob;

Gpt        = N * pts;                  % GP in global coord, used


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
        nodeId = sctr(in);
        enrnoI = enrich_node(nodeId); 
        dRidx  = dRdx(in,1);
        dRidy  = dRdx(in,2);
        Ri     = N(in);
        
        if ( enrnoI == 1)     % H(x) enriched node
            
            crackId = crack_node(nodeId);
            
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            
            % Enrichment function, H(x) at global Gauss point
            dist = signed_distance(xCr,Gpt);
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            dist = signed_distance(xCr,controlPts(nodeId,:));
            Hi   = heaviside(dist);
            %Hi=0;
            % Bxfem at node "in"
            
            aa   = dRidx*(Hgp - Hi);
            bb   = dRidy*(Hgp - Hi);
            
            BI_enr = [aa 0 ;
                      0 bb;
                      bb aa];
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
        elseif ( enrnoI == 2) % B(x) enriched node
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = branch(r,theta,alpha);
            
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta); %BrI = zeros(4);
            
            % components of Benr matrix
            
            aa = dRidx*(Br(1)-BrI(1)) + Ri*dBdx(1) ;
            bb = dRidy*(Br(1)-BrI(1)) + Ri*dBdy(1) ;
            B1_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(2)-BrI(2)) + Ri*dBdx(2) ;
            bb = dRidy*(Br(2)-BrI(2)) + Ri*dBdy(2) ;
            B2_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(3)-BrI(3)) + Ri*dBdx(3) ;
            bb = dRidy*(Br(3)-BrI(3)) + Ri*dBdy(3) ;
            B3_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(4)-BrI(4)) + Ri*dBdx(4) ;
            bb = dRidy*(Br(4)-BrI(4)) + Ri*dBdy(4) ;
            B4_enr = [aa 0 ; 0 bb ; bb aa];
            
            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
        elseif ( enrnoI == 3)     % material interface enriched node
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                %E    = abs(Zm)-abs(chi(in)); % shift enrichment
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; % derivatives w.r.t global coords
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            
        elseif ( enrnoI == 4) % bi-mat crack tip/mat interface enriched node
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = biMatCrackBranch(r,theta,ep,alpha);
                      
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = biMatCrackBranchNode(r,theta,ep,alpha);
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa = dRidx*(Br(i)) + Ri*dBdx(i) ;
                bb = dRidy*(Br(i)) + Ri*dBdy(i) ;
                BI_enr = [aa 0 ; 0 bb ; bb aa];
                                
                Bxfem = [Bxfem BI_enr];                                
            end
            
            % enrichment functions and derivatives due to material interface 
            
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                %E    = abs(Zm)-abs(chi(in)); % shift enrichment
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; 
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            Bxfem = [Bxfem BI_enr]; 
            clear BI_enr ;

         elseif ( enrnoI == 5) % bi-mat crack tip enriched node
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = biMatCrackBranch(r,theta,ep,alpha);
                      
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = biMatCrackBranchNode(r,theta,ep,alpha);
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa = dRidx*(Br(i)) + Ri*dBdx(i) ;
                bb = dRidy*(Br(i)) + Ri*dBdy(i) ;
                BI_enr = [aa 0 ; 0 bb ; bb aa];
                                
                Bxfem = [Bxfem BI_enr];                                
            end
        elseif ( enrnoI == 6) % 4 B(x) enriched node and 1 weak enr. func.
            crackId = crack_node(nodeId);
            xTip    = xTips(crackId,:);
            xCr     = reshape(xCrack(crackId,:,:),2,2);
            seg     = xCr(2,:) - xCr(1,:);   % tip segment
            alpha   = atan2(seg(2),seg(1));  % inclination angle
            QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [Br,dBdx,dBdy] = branch(r,theta,alpha);
            
            % compute branch functions at node "in"
            
            xp    = QT*(controlPts(nodeId,:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            
            if ( theta > pi || theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            
            [BrI] = branch_node(r,theta); %BrI = zeros(4);
            
            % components of Benr matrix
            
            aa = dRidx*(Br(1)-BrI(1)) + Ri*dBdx(1) ;
            bb = dRidy*(Br(1)-BrI(1)) + Ri*dBdy(1) ;
            B1_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(2)-BrI(2)) + Ri*dBdx(2) ;
            bb = dRidy*(Br(2)-BrI(2)) + Ri*dBdy(2) ;
            B2_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(3)-BrI(3)) + Ri*dBdx(3) ;
            bb = dRidy*(Br(3)-BrI(3)) + Ri*dBdy(3) ;
            B3_enr = [aa 0 ; 0 bb ; bb aa];
            
            aa = dRidx*(Br(4)-BrI(4)) + Ri*dBdx(4) ;
            bb = dRidy*(Br(4)-BrI(4)) + Ri*dBdy(4) ;
            B4_enr = [aa 0 ; 0 bb ; bb aa];
            
            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
            
            % enrichment functions and derivatives due to material interface 
            
            chi  = CHI(sctr);                        
            Zm   = dot(N,chi);                        
            
            % enrichment functions and derivatives
            
            if     (weakEnrFunc == 1)
                % Moes function
                achi = abs(chi);
                Za   = dot(N,achi);
                Zmn  = Zm/abs(Zm);
                E    = Za-abs(Zm);
                Exi  = dot(dRdxi,achi)  - Zmn*dot(dRdxi,chi);
                Eet  = dot(dRdeta,achi) - Zmn*dot(dRdeta,chi);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
                Exi  = sign(Zm)*dot(dRdxi,chi);
                Eet  = sign(Zm)*dot(dRdeta,chi);
            end
            
            aa = dRdxi (in)*E + Ri*Exi;
            bb = dRdeta(in)*E + Ri*Eet;
            
            t  = [aa bb]*invJacob; 
            
            BI_enr = [t(1) 0;0 t(2);t(2) t(1)] ;
            
            Bxfem = [Bxfem BI_enr]; 
            clear BI_enr ;
        end
    end          % end of loop on nodes
    % B matrix
    B = [Bfem Bxfem];
    clear Bfem; clear Bxfem;
end              % end of switch between enriched and non-enriched elements

