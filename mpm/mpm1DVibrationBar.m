% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with material points.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized
% Material Point Method", Buzzi et al, CMES 2008.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.

%%

addpath ../fem_util/


%%
clc
clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%%
%
E = 100;               % Young modulus
L = 25;                % length of the bar
rho = 1;               % density

%%
%  Computational grid: two-noded elements

elemCount = 63;
nodes     = linspace(0,L,elemCount+1);
elements  = zeros(elemCount,2);

elemCount = elemCount + 1;

for ie=1:elemCount
    elements(ie,:) = ie:ie+1;
end

nodes = [nodes L+1.5];
elements = [elements;[14 15]];
nodeCount = size(nodes,2);

velo = zeros(nodeCount,1);
mId  = floor((elemCount-1)/2)+1;

%%
% Material points: centers of the elements
% one particle per element

c      = sqrt(E/rho);
beta1  = pi/2/L;
omega1 = beta1*c;
deltax = L/elemCount;

xp = zeros(elemCount-1,1);

for p=1:elemCount-1
    xp(p) = 0.5*(nodes(p) + nodes(p+1));
end

pCount = length(xp);

Mp = L/(elemCount-1)*ones(pCount,1);     % mass
Vp = L/(elemCount-1)*ones(pCount,1);     % volume
Fp =             ones(pCount,1);     % gradient deformation
Vp0=Vp;

s   = zeros(pCount,1);               % stress
eps = zeros(pCount,1);               % strain
vp  = zeros(pCount,1);               % velocity

% initial velocities

for p=1:pCount
    vp(p) = 0.1*sin(beta1*xp(p));
end

%q  = Mp*vp;                          % momentum

%%
hold on
plot(nodes,zeros(nodeCount,1)+1/2,'r-s');
plot(xp,zeros(pCount,1)+1/2,'b*');
axis([0 L+0.5 0 1.])

%%
% data structure to store the material points for each element
% this data structure is updated for every time step

for e=1:elemCount
    mpoints{e}=e;
end
mpoints{elemCount}=[];

%% Time loop

tol = 1e-5;


dtime = 0.1*deltax/c;
time  = 100;
t     = 0;

ta = [];           % time
va = [];           % velocities 
xa = [];           % position
ka = [];           % kinetic energy 
sa = [];           % strain energy

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,1);  % nodal momentum vector
niforce   = zeros(nodeCount,1);  % nodal internal force vector
neforce   = zeros(nodeCount,1);  % nodal external force vector

while ( t < time )
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        Le    = enode(2)-enode(1);
        mpts  = mpoints{e};
        
        % loop over particles
        
        for p=1:length(mpts)
            pid  = mpts(p);
            
            N1  = 1 - abs(xp(pid)-enode(1))/Le;
            N2  = 1 - abs(xp(pid)-enode(2))/Le;
            
            dN1 = -1/Le;
            dN2 =  1/Le;
            
            % particle mass and momentum to node
            
            nmass(esctr(1))      = nmass(esctr(1))     + N1*Mp(pid);
            nmass(esctr(2))      = nmass(esctr(2))     + N2*Mp(pid);
            nmomentum(esctr(1))  = nmomentum(esctr(1)) + N1*Mp(pid)*vp(pid);
            nmomentum(esctr(2))  = nmomentum(esctr(2)) + N2*Mp(pid)*vp(pid);
            
            % internal force
            
            niforce(esctr(1)) = niforce(esctr(1)) - Vp(pid)*s(pid)*dN1;
            niforce(esctr(2)) = niforce(esctr(2)) - Vp(pid)*s(pid)*dN2;
        end
    end
    
    % update nodal momenta
    
    nmomentum(1)  = 0; % Boundary conditions f1 = m1*a1, a1=0
    niforce(1)    = 0;
    
    for i=1:nodeCount
        nmomentum(i) = nmomentum(i) + niforce(i)*dtime;
    end
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = elements(e,:);
        enode = nodes(esctr);
        Le    = enode(2)-enode(1);
        mpts  = mpoints{e};
        
        % loop over particles
        
        for p=1:length(mpts)
            pid  = mpts(p);
            
            N1  = 1 - abs(xp(pid)-enode(1))/Le;
            N2  = 1 - abs(xp(pid)-enode(2))/Le;
            
            dN1 = -1/Le;
            dN2 =  1/Le;
            
            if nmass(esctr(1)) > tol 
               vp(pid)  = vp(pid) + dtime * N1*niforce(esctr(1))/nmass(esctr(1));                
            end
            
            if nmass(esctr(2)) > tol 
               vp(pid)  = vp(pid) + dtime * N2*niforce(esctr(2))/nmass(esctr(2));                
            end
            
            xp(pid)  = xp(pid) + dtime * (N1*nmomentum(esctr(1))/nmass(esctr(1)) + ...
                 N2*nmomentum(esctr(2))/nmass(esctr(2)));

            v1 = nmomentum(esctr(1))/nmass(esctr(1));
            v2 = nmomentum(esctr(2))/nmass(esctr(2));
    
            %if ( esctr(1) == 1 ) v1 = 0; end
            % gradient velocity
            
            Lp = dN1 * v1 + dN2 * v2;
            
            Fp(pid) = (1 + Lp*dtime)*Fp(pid);
            Vp(pid) = Fp(pid)*Vp0(pid);                                    
            dEps    = dtime * Lp;                                    
            s(pid)  = s(pid)   + E * dEps;    
            eps(pid)= eps(pid) + dEps;  
            
            k = k + 0.5*vp(pid)^2*Mp(pid);
            u = u + 0.5*s(pid)*eps(pid)*Mp(pid);
        end
    end
%     velo(:) = 0;
%     for e=1:elemCount
%         esctr = elements(e,:);
%         enode = nodes(esctr);
%         Le    = enode(2)-enode(1);
%         mpts  = mpoints{e};
%         
%         % loop over particles
%         
%         for p=1:length(mpts)
%             pid  = mpts(p);
%             
%             N1  = 1 - abs(xp(pid)-enode(1))/Le;
%             N2  = 1 - abs(xp(pid)-enode(2))/Le;
%             
%             dN1 = -1/Le;
%             dN2 =  1/Le;
%             
%             velo(esctr(1)) = velo(esctr(1)) + N1*Mp(pid)*vp(pid)/nmass(esctr(1));
%             velo(esctr(2)) = velo(esctr(2)) + N2*Mp(pid)*vp(pid)/nmass(esctr(2));
%             %         if (nmass(esctr(1)) > tol ) v1 = nmomentum(esctr(1))/nmass(esctr(1));
%             %         else v1 = 0; end
%             %         if (nmass(esctr(2)) > tol ) v2 = nmomentum(esctr(2))/nmass(esctr(2));
%             %         else v2 = 0; end
%         end
%     end
%     for e=1:elemCount
%         esctr = elements(e,:);
%         enode = nodes(esctr);
%         Le    = enode(2)-enode(1);
%         
%         mpts  = mpoints{e};
%         
%         % loop over particles
%         
%         for p=1:length(mpts)
%             pid  = mpts(p);
%             
%             N1  = 1 - abs(xp(pid)-enode(1))/Le;
%             N2  = 1 - abs(xp(pid)-enode(2))/Le;
%             
%             dN1 = -1/Le;
%             dN2 =  1/Le;
%             velo(1)=0;
%             % gradient velocity
%             
%             Lp = dN1 * velo(esctr(1)) + dN2 * velo(esctr(2));
%             
%             Fp(pid) = (1 + Lp*dtime)*Fp(pid);
%             Vp(pid) = Fp(pid)*Vp0(pid);
%             
%             % strain increment
%             
%             dEps = dtime * Lp;
%             
%             % stress update
%             
%             s(pid) = s(pid) + E * dEps;
%             
%             if nmass(esctr(1)) > tol
%                 xp(pid)  = xp(pid) + dtime * N1*nmomentum(esctr(1))/nmass(esctr(1));
%             end
%             
%             if nmass(esctr(2)) > tol
%                 xp(pid)  = xp(pid) + dtime * N2*nmomentum(esctr(2))/nmass(esctr(2));
%             end
%             
%         end
%     end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    va = [va;vp(mId)];
    xa = [xa;xp(mId)];    
    ka = [ka;k];
    sa = [sa;u];
    
    % update the element particle list
    
    
    pe = round(xp/deltax);
    
    
    for e=1:elemCount
        id  = find(pe==e);
        %mpoints{e}=id;
    end
    
    % advance to the next time step
    
    t = t + dtime;
end

%%

% exact solution

vExact = 0.1*cos(omega1.*ta)*sin(beta1*0.5*L);
%vExact = 0.1/(beta1*L)*cos(omega1.*ta);
uExact = (0.1/omega1)*sin(omega1.*ta)*sin(beta1*0.5*L);

% % rho    = 1;
% w      = 1/L*sqrt(E/rho);
% vExact = 0.1*cos(w.*ta);
% xExact = 0.5*exp((0.1/L/w)*sin(w.*ta));

figure
set(gca,'FontSize',14)
hold on
plot(ta,va,'b-','LineWidth',1.6);
plot(ta,vExact,'r--','LineWidth',2);
xlabel('Time')
ylabel('Velocity')
legend('MPM','Exact')
axis([0 100 -0.1 0.1])

figure
set(gca,'FontSize',14)
hold on
plot(ta,xa-0.5*L,'b-','LineWidth',1.6);
plot(ta,uExact,'r--','LineWidth',2);
xlabel('Time')
ylabel('Displacement')
legend('MPM','Exact')
set(gca,'FontSize',14)
axis([0 100 -0.15 0.2])

figure 
set(gca,'FontSize',14)
hold on
plot(ta,ka,'b-','LineWidth',1.6);
plot(ta,sa,'r--','LineWidth',2);
plot(ta,ka+sa,'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 100 0 0.08])

