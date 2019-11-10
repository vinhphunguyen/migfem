% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotStress(stressXX,stressXY,stressYY,element,node)
% This function plots the stress for quadrilateral elements calculated in
% elemStress.m.

nElem = size(element,1);

for iPlot = 1:3
    if iPlot == 1
        Es = stressXX;                                                      % Plot sigmaXX
    elseif iPlot == 2
        Es = stressXY;                                                      % Plot sigmaXY
    elseif iPlot == 3
        Es = stressYY;                                                      % Plot sigmaYY 
    end
    
    figure; hold on;
    
    smin = min(min(Es));                                                    % Find the minimum stress value
    smax = max(max(Es));                                                    % Find the maximum stress value

    if (smin < 1E-3) && (smax < 1E-3)                                       % Check if stresses are simply numerical noise
        if size(Es,2) == 1                                                  % Averaged stress values
            Es = zeros(length(Es));                                         % Define stress to be zero
        elseif size(Es,2) == 4                                              % Actual stress values
            Es = zeros(length(Es),4);                                       % Define stress to be zero
        end
        caxis([-0.1 0.1]);                                                  % Define bounds if stress is zero
    elseif (smin - smax) < 1E-6
        caxis([smin-0.1*smin smax+0.1*smax]);                               % Define bounds if stress is nonzero, constant
    else
        caxis([smin smax]);                                                 % Define bounds if stress is not constant
    end

    axis equal; axis off; colorbar('horiz');
    
    Ex = zeros(4,nElem); Ey = Ex; Ec = Ex;                                  % Initialize plotting matrices

    for iElem = 1:nElem
        NN          = element(iElem,:);                                    % Nodes of current element
        Ex(:,iElem) = node(NN',1);                                           % X-coordinates of nodes
        Ey(:,iElem) = node(NN',2);                                           % Y-coordinates of nodes
        if size(Es,2) == 1                                                  % Averaged stress values
            Ec(:,iElem) = Es(NN');                                          % Stress values at nodes
        elseif size(Es,2) == 4                                              % Actual stress values
            Ec(:,iElem) = Es(iElem,:)';                                     % Stress values at nodes
        end
    end
    
    % patch(Ex,Ey,Ec);
    H = patch(Ex,Ey,Ec);                                                    % Plot the stresses
    set(H,'LineStyle','none')
     
    if iPlot == 1
        title('\sigma_X_X');
    elseif iPlot == 2
        title('\sigma_X_Y');
    elseif iPlot == 3
        title('\sigma_Y_Y');
    elseif iPlot == 4
        title('\sigma_V_M');
    end
           
    hold off
end