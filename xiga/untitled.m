% =======================================================================
% The original version is written by Stephane Bordas
% Modified by Nguyen Vinh Phu 2005.05.26 for XFEM- Homogenous
%             Truong Quang Tri 2006 for XFEM - interfacial crack
% EMMC IX, Hochiminh University of Technology, Vietnam
% ======================================================================= 

clear all; close all; clc;
format long;

addpath('~/code/xfem-efg-matlab/fem_util');

% Data file for edge cracked plate in tension
% using uniform Q4 mesh

% Geometry parameters

L = 3 ;
W = 9;
% Corner points of the rectangle domain

pt1 = [ 0 0];
pt2 = [ L 0];
pt3 = [ L W];
pt4 = [ 0 W];

nnx = 60;
nny = 180;

[node,element] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,'Q4'); 

% define crack

xcr=[-0.01 0.5*L];
ycr=[1.5*L 1.5*L];

%define line interface
xc = 0.0;
yc = 0.0;
xi = [-0.01 L];
yi = [1.5*L 1.5*L];

% DEFINE BOUNDARIES
% Here we define the boundary discretizations.
uln=nnx*(nny-1)+1;
urn=nnx*nny;
lrn=nnx;
lln=1;                          % lower left node number
cln=nnx*(nny-1)/2+1;  % node number at (0,0)

leftEdge  = [ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
rightEdge = [ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
topEdge   = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge   = [ lln:1:(lrn-1); (lln+1):1:lrn ]';

leftNodes   = unique(leftEdge);
rightNodes  = unique(rightEdge);
topNodes    = unique(topEdge);
botNodes    = unique(botEdge);

dispx      = botNodes;
dispy      = botNodes;
dispNodes  = unique([dispx; dispy]);
tractNodes = unique(topNodes) ;

nodalLoad  = L/(nnx-1);

% compute the size of domain
numnode = size(node,1) ;
numelem = size(element,1) ;

% Recognize elt type automatically
switch size(element,2)
 case 3
  elt_type = 'T3';
 case 4
  elt_type = 'Q4';
 otherwise
  error('Unknown element type found in mesh');
end

% ================================================================
% plot the mesh to check the BCs are correct
colordef black
plot_mesh(node,element,elt_type,'g.-');
hold on
plot(node(dispNodes,1),node(dispNodes,2),'r>');
plot(node(tractNodes,1),node(tractNodes,2),'r*');

inter=plot(xi,yi,'b');
cr = plot(xcr,ycr,'y.');
set(cr,'LineWidth',4);
axis off
pause


% *************************************************************************
% MATERIAL CONSTANT
% *************************************************************************

E1  = 200.0 ; % Young modulus
nu1 = 0.3 ; % Poisson ratio

E2  = 100.0; %Youngmodulus
nu2 = 0.3 ; % Poisson ratio

% in the following, 1, 2, and 3 refer to the 3 regions

elemType1 = elt_type ; % get the element type from the Gmsh file

% BOUNDARY CONDITIONS
dispval = 0.0 ; % value of fixed displacement

% LOADING
delta = 0.0; % prescribe displacement

% -------------------------------------------------------------
% WRITE THE INPUT FILE follows format of FEMOBJ
% For more info on FEMOBJ, see www.zace.com
% -------------------------------------------------------------
outFileName = '~/Downloads/openxfem++-1/temp.m';
outfid = fopen(outFileName,'w');
fprintf(outfid,'* Solving the Griffith problem using the XFEM *\n');
fprintf(outfid,'* The XFEM C++ library developped by *\n');
fprintf(outfid,'* S. Bordas C. Dunant N.V.Phu T.Q.Tri R. Duddu *\n');
fprintf(outfid,'* ---------------------------------------------------------\n'); 
% Export result to matlab file
fprintf(outfid,'\nExportDisplacement 1 *\n');
fprintf(outfid,'ExportStress 1 *\n');
fprintf(outfid,'ExportStrain 1 *\n\n');

% -----------------------------------
% MATERIAL
% -----------------------------------

fprintf(outfid,'Material 2 \n');
fprintf(outfid,'1 class ElasticMaterial E %15.12f n %15.12f t 1.0 *\n',E1,nu1); 
fprintf(outfid,'2 class ElasticMaterial E %15.12f n %15.12f t 1.0 *\n',E2,nu2); 
fprintf(outfid,'**\n\n');

% ---------------------------------------
% TIME INTEGRATION
% ---------------------------------------

fprintf(outfid,'TimeIntegrationScheme *\n');
fprintf(outfid,'1 class Static *\n');
fprintf(outfid,'**\n\n');

% ------------------------------------------------------------

% LOAD SECTION (includes Dirichlet conditions)
% Example :
% Load 2
% 1 class BoundaryCondition loadTimeFunction 1 conditions 1 d 0. *
% 2 class NodalLoad loadTimeFunction 2 components 2 4000. 0. *
% **
% ------------------------------------------------------------

fprintf(outfid,'Load %d \n',3);
fprintf(outfid,'%d class BoundaryCondition loadTimeFunction 1 conditions 1 d %15.14f  * \n',...
    1,0);
fprintf(outfid,'%d class NodalLoad   loadTimeFunction 1 components  2  %16.14f  %16.14f  * \n',...
    2,0,nodalLoad );
fprintf(outfid,'%d class NodalLoad   loadTimeFunction 1 components  2  %16.14f  %16.14f  * \n',...
    3,0,0.5*nodalLoad );

fprintf(outfid,'**\n\n');

% -------------------------------------------------------------------------
%           LOADTIME FUNCTIONS
% -------------------------------------------------------------------------

fprintf(outfid,'LoadTimeFunction 1 \n');
fprintf(outfid,'1 class ConstantFunction f(t) %f *\n',1);
fprintf(outfid,'**\n\n');

% -------------------------------------------------------------------------------
% NODE SECTION
% -------------------------------------------------------------------------------

fprintf(outfid,'Node %d \n',numnode);

% FIRST WRITE THE FINITE ELEMENT NODES ASSOC WITH ELEMENTGS IN DOMAIN 1
nodei = 0;

% THEN WRITE THE FINITE ELEMENT NODES ASSOC WITH ELEMENTGS IN DOMAIN 3
for i = 1 : numnode
    nodei = nodei +1 ;
    % always written
    fprintf(outfid,['%d nDofs 2 coord 2 %16.12f %16.12f'],nodei,node(i,1),node(i,2));
    
    if (ismember(i,dispx))
        fprintf(outfid,' bcOnDof1 1 ');     % clamped nodes
    end
    
    if (ismember(i,dispy))
        fprintf(outfid,' bcOnDof2 1 ');
    end
    
    if (ismember(i,tractNodes))
        fprintf(outfid,' loads 1 2 ');
    end
    fprintf(outfid,' *');
    fprintf(outfid,' \n');
end

fprintf(outfid,'**\n\n');

% -------------------------------------------------------------------------
% % % % ELEMENT SECTION
%
%
%
% ------------------------------------------------------------------------- 
% % COMPUTE THE TOTAL NUMBER OF ELEMENTS

fprintf(outfid,'Element %d \n',numelem);
nume = 0;
% EXPORT FINTE ELEMENTS SECOND (THIS IS MANDATORY)

for i=1:numelem
    nume = nume+1;
    fprintf(outfid,'%d class Q4U nodes %d %d %d %d  mat 1 ',nume,...
        element(i,1),element(i,2),element(i,3),element(i,4));
    fprintf(outfid,' *');
    fprintf(outfid,' \n');
end
fprintf(outfid,'**\n\n');

% ----------------------------------------------------------------------------------
% % % % ENRICHMENT ITEM SECTION
% Example :  for one crack
% EnrichmentItem 3
% 1 class CrackInterior geometry 1  EnrichmentFunctions 1 1  enrichScheme 3*
% 2 class CrackTip Type BimatElast Mat 2 1 2 geometry 2  EnrichmentFuncs  4 2  3  4  5  enrichScheme 1 enrichRadius 2.0 *
% 3 class CrackTip Type BimatElast Mat 2 1 2 geometry 3  EnrichmentFuncs  4 2 3 4 5 enrichScheme 1 enrichRadius 2.0 *
% **
% -------------------------------------------------------------------------- 

numOfEnrItem = 3;
enrichScheme = 1 ; % standard enrichment => no need enrichRadius
enrichRadius = 3 ; % radius of enrichment
domainIntFact = 2.5 ; % the radius of the domain will be domainIntRadius * sqrt(area_of_tipElement)

disp('ATTENTION! JUST FOR ONE CRACK !!! TOO BAD');

fprintf(outfid,'EnrichmentItem %d \n',numOfEnrItem);
fprintf(outfid,'%d class CrackInterior  geometry 1 myTips 1 2 EnrichmentFuncs 1  1  enrichScheme 3',1);
fprintf(outfid,' *');
fprintf(outfid,' \n');

if (enrichScheme == 1)
   fprintf(outfid,'%d class CrackTip  Type BiMatElast Mat 2  1 2 geometry %d EnrichmentFuncs  12  2 3 4 5 6 7 8 9 10 11 12 13 enrichScheme %d domainIntFac 5.3f ',...
        2,3,enrichScheme,domainIntFact); 
   fprintf(outfid,' *'); 
   fprintf(outfid,' \n');
else
    fprintf(outfid,'%d class CrackTip  Type BiMatElast Mat 2  1 2  geometry %d EnrichmentFuncs  12  2 3 4 5 6 7 8 9 10 11 12 13 enrichScheme %d domainIntFac 5.3f domainIntRadius %5.3f ',...
        3,nsegc+2,enrichScheme,enrichRadius,domainIntFact);
    fprintf(outfid,' *'); 
    fprintf(outfid,' \n');
end

fprintf(outfid,'%d class MaterialInterface  geometry %d Mat 2 1 2 EnrichmentFuncs 1 14  enrichScheme 3',...
    3,4);
fprintf(outfid,' *');
fprintf(outfid,' \n');

fprintf(outfid,'**\n\n');

% -----------------------------------------------------------------------
% % GEOMETRY ENTITY SECTION
% Example :
%  GeometryEntity 3
%  1 class PiecewiseLinear numOfVertices 2 vertices 2   3
%  2 class Vertex coord 0.75 3.0                    *
%  3 class Vertex coord 1.25 3.0                    *
%  **
% ----------------------------------------------------------------------- 

numOfGeoEntity = 6;
fprintf(outfid,'GeometryEntity %d \n',numOfGeoEntity);

%define crack
fprintf(outfid,'%d class PiecewiseLinear  numOfVertices %d vertices ',1,2);
for i=1:2
    fprintf(outfid,'%d ',i+1);
end
fprintf(outfid,'geoDescription 1 *\n');

for i=1:2
    fprintf(outfid,'%d class Vertex  coord 2 %12.9f %12.9f geoDescription 1',...
                    i+1,xcr(i),ycr(i));
    fprintf(outfid,' *');fprintf(outfid,' \n');
end

%define material interface
fprintf(outfid,'%d class PiecewiseLinear  numOfVertices 2 vertices 5 6 geoDescription 1',...
                4); 
fprintf(outfid,' *');fprintf(outfid,' \n');
fprintf(outfid,'%d class Vertex  coord 2  %4.2f %4.2f geoDescription 1',...
                  5,xi(1),yi(1));  
fprintf(outfid,' *');fprintf(outfid,' \n');
fprintf(outfid,'%d class Vertex  coord 2  %4.2f %4.2f geoDescription 1',...
                  6,xi(2),yi(2));                
fprintf(outfid,' *'); 
fprintf(outfid,' \n');
fprintf(outfid,'**\n\n');

% --------------------------------------------------------
%
% ENRICHMENT FUNCTION SECTION
%
% --------------------------------------------------------

numOfEnrFunc = 14;
fprintf(outfid,'EnrichmentFunction %d \n',numOfEnrFunc);
fprintf(outfid,'%d class DiscontinuousField ',1); 
fprintf(outfid,' *'); 
fprintf(outfid,' \n');

for i = 1:12
    fprintf(outfid,['%d class BiMatCrackAsymp' num2str(i) ],1+i); 
    fprintf(outfid,' *'); 
    fprintf(outfid,' \n');
end

fprintf(outfid,'%d class AbsSignedDistance ',14); 
fprintf(outfid,' *'); 
fprintf(outfid,' \n');
fprintf(outfid,'**\n\n');

% -----------------------------------
% NLSOLVER
% -----------------------------------

fprintf(outfid,'NLSolver *\n');
fprintf(outfid,'1 class NewtonRaphson n 100 t 1e-5 c 1 *\n');
fprintf(outfid,'**\n\n');

% -----------------------------------
% TIME STEP
% -----------------------------------

fprintf(outfid,'TimeStep 1 \n');
fprintf(outfid,'1 dt 1.0 *\n');
fprintf(outfid,'**\n\n');

% -----------------------------------
% CLOSE FILE
% -----------------------------------

fclose(outfid);

disp('---------------------------------------------------------------');
disp('number of elements in set ');
length(element)
disp('number of nodes');
length(node)
disp('---------------------------------------------------------------');

% END OF SCRIPT
% ========================================================================

plot_field(node+fac*[U_1(:,1) U_1(:,2)],element,'Q4',U_1(:,1));
colorbar

plot_field(node+fac*[U_1(:,1) U_1(:,2)],element,'Q4',sigma_(:,1));
colorbar