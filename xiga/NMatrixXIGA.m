function [exN] = NMatrixXIGA(e,enrich_node,N)
% compute extended shape functions
% exN = [N1 N2 N3 | N1(H(x)-H1) N2(H(x)-H2) N3(phi(x)-phi3)]
% for elements with three nodes, two of them are enriched
% with Heaviside function and one is enriched with phi(x).
% Usage:
% This function is used to determine the displacement
% at a point i.e. in displacement visualization.
% Vinh Phu Nguyen
% Johns Hopkins University

global controlPts element xCrack xTips crack_node CHI ep weakEnrFunc

sctr = element(e,:);
nn   = length(sctr);
pts  = controlPts(sctr,:);
Gpt  = N * pts;                  % GP in global coord, used

% Nfem is always computed

Nfem = N;

% Switch between non-enriched and enriched elements
if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    exN = Nfem;
else                               % Enriched elements
    Nxfem = [] ;
    % loop on nodes, check node is enriched ...
    for in = 1 : nn
        nodeId = sctr(in);
        enrnoI = enrich_node(nodeId);
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
            % Nxfem at node "in"
            
            aa   = Ri*(Hgp - Hi);
            
            % Add to the total Nxfem
            Nxfem = [Nxfem aa];
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
            
            [BrI] = branch_node(r,theta);%BrI = zeros(4);
            
            % components of Nxfem vector
            
            aa1 = Ri*(Br(1)-BrI(1));
            aa2 = Ri*(Br(2)-BrI(2));
            aa3 = Ri*(Br(3)-BrI(3));
            aa4 = Ri*(Br(4)-BrI(4));
            
            Nxfem = [Nxfem aa1 aa2 aa3 aa4];
        elseif ( enrnoI == 3) % inclusion enriched node
            chi  = CHI(sctr);
            achi = abs(chi);
            Zm   = dot(N,chi);
            Za   = dot(N,achi);
            
            if     (weakEnrFunc == 1)
                E    = Za-abs(Zm);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
            end
           
            % Add to the total Nxfem
            aa   = Ri*E;
            Nxfem = [Nxfem aa];
            
        elseif ( enrnoI == 4) % bi-mat crack tip enriched node
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
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa    = Ri*Br(i);                
                Nxfem = [Nxfem aa];
            end
            
            % due to material interface
            
            chi  = CHI(sctr);
            achi = abs(chi);
            
            Zm   = dot(N,chi);
            Za   = dot(N,achi);
            
            % enrichment functions
            
            if     (weakEnrFunc == 1)
                E    = Za-abs(Zm);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
            end
            
            aa = Ri*E;
            Nxfem = [Nxfem aa];
        
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
            
            % components of Benr matrix due to branch functions
            
            for i=1:12
                aa    = Ri*Br(i);                
                Nxfem = [Nxfem aa];
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
            
            [BrI] = branch_node(r,theta);%BrI = zeros(4);
            
            % components of Nxfem vector
            
            aa1 = Ri*(Br(1)-BrI(1));
            aa2 = Ri*(Br(2)-BrI(2));
            aa3 = Ri*(Br(3)-BrI(3));
            aa4 = Ri*(Br(4)-BrI(4));
            
            Nxfem = [Nxfem aa1 aa2 aa3 aa4];
            
            % due to material interface
            
            chi  = CHI(sctr);
            achi = abs(chi);
            
            Zm   = dot(N,chi);
            Za   = dot(N,achi);
            
            % enrichment functions
            
            if     (weakEnrFunc == 1)
                E    = Za-abs(Zm);
            elseif (weakEnrFunc == 2)
                E    = abs(Zm);
            end
            
            aa = Ri*E;
            Nxfem = [Nxfem aa];
        end
    end          % end of loop on nodes
    
    % N vector
    exN = [Nfem Nxfem];
    clear Nxfem;
    clear Nfem;
end              % end of switch between enriched and non-enriched elements

