function [K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed)

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';
f=f-K(:,wdofs)*wFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(wdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(:,wdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(wdofs,wdofs)=bcwt*speye(length(wdofs));
