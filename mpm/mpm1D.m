% This file implements the Material Point Method of Sulsky 1994.
% One dimensional problem with one material point.
% The grid is one two-noded linear element.
% This example is taken from "Caveats on the Implementation of the Generalized 
% Material Point Method", Buzzi et al, CMES 2008.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.

%%
addpath ../fem_util/

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

%%
clc
clear all;

%%
%
E = 4*pi^2;            % Young modulus
L = 1;                 % length of the bar
rho    = 1;
%%
%  Computational grid

nodes     = [0 L];
elements  = [1 2];

elemCount = size(elements,1);
nodeCount = size(nodes,2);

%%
% Material points

xp = 0.5*L;            % position
Mp = 1;                % mass 
Vp = 1;                % volume
vp = 0.1;              % velocity
s  = 0.;               % stress
eps=0;                 % strain 
q  = Mp*vp;            % momentum

%% Time loop

c     = sqrt(E/rho);
cfl   = c/L;
time  = 2;
dtime = 0.001;
t     = 0;

ta = [];    % time
va = [];    % velocities
xa = [];    % positions
ka = [];    % kinetic energy
sa = [];    % strain energy

while ( t < time )
    
    % shape functions and derivatives
    
    N1  = 1 - abs(xp-nodes(1))/L;
    N2  = 1 - abs(xp-nodes(2))/L;
    
    dN1 = -1/L;
    dN2 =  1/L;
    
    % particle mass and momentum to node
    
    m1   = N1*Mp;
    mv1  = N1*q;
    
    m2   = N2*Mp;
    mv2  = N2*q;
    
    mv1  = 0; % Boundary condition v=0 at node 1
    
    % internal force
    
    fint1 = - Vp*s*dN1;
    fint2 = - Vp*s*dN2;
    
    f1    = fint1 ;   % there is no external force
    f2    = fint2 ;
    
    % update nodal momenta

    f1  = 0; % Boundary conditions f1 = m1*a1, a1=0
    
    mv1 = mv1 + f1*dtime;
    mv2 = mv2 + f2*dtime;
    
    % update particle velocity and position
    
    vp  = vp + dtime * (N1*f1/m1 + N2*f2/m2);
    xp  = xp + dtime * (N1*mv1/m1 + N2*mv2/m2);
    
    q   = Mp*vp;
    
    v1 = N1*Mp*vp/m1;
    v2 = N2*Mp*vp/m2;
    v1 = 0;
    % gradient velocity
    
    Lp = dN1 * v1 + dN2 * v2;
    
    % strain increment
    
    dEps = dtime * Lp;
    
    % stress update
    
    eps  = eps + dEps;
    s    = s + E * dEps;
            
    % store time,velocty for plotting
    
    ta = [ta;t];
    va = [va;vp];
    xa = [xa;xp];
    
    ka = [ka;0.5*vp*vp*Mp];
    sa = [sa;0.5*s*eps*Mp];
    
    % advance to the next time step
    
    t = t + dtime;
end

%%

% exact solution


w      = 1/L*sqrt(E/rho);
vExact = 0.1*cos(w.*ta);
xExact = 0.5*exp((0.1/L/w)*sin(w.*ta));

figure 
set(gca,'FontSize',14)
hold on
plot(ta,va,'b-','LineWidth',1.6);
plot(ta,vExact,'r--','LineWidth',2);
%set(gca,'YTick',va)
%set(gca,'YTickLabel',sprintf('%4.2f|',va))
%set(gca,'YTick',[-0.15 -0.10 -0.05 0.00 0.05 0.10 0.15])
xlabel('Time')
ylabel('Velocity')
legend('MPM','Exact')


figure 
set(gca,'FontSize',14)
hold on
plot(ta,xa-0.5,'b-','LineWidth',1.6);
plot(ta,xExact-0.5,'r--','LineWidth',2);
xlabel('Time')
ylabel('Displacement')
legend('MPM','Exact')

figure 
set(gca,'FontSize',14)
hold on
plot(ta,ka,'b-','LineWidth',1.6);
plot(ta,sa,'r--','LineWidth',2);
plot(ta,ka+sa,'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
axis([0 2 0 0.007])

