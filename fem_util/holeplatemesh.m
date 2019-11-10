function a=holeplatemesh(cpt,w,h,r,oct,nr,nt,ratior,ratiot)
  
%  HOLEPLATEMESH(center pt,w,h,r,octant,nr,nt,ratior,ratiot)
%
%  Creates an ocatant mesh of a hole in a plate
  

switch(oct)
 case 1
  
  apts=linemesh([r,0],[w/2,0],nr,ratior);
  
  theta1=0;
  theta2=atan2(h,w);
   
  bpts=linemesh([cos(theta2),sin(theta2)]*r,[w/2,h/2],nr,ratior);
  
 case 2
   
 case 3
   
 case 4 

 case 5
  
 case 6
  
 case 7
 
 case 8
  
end

theta=linemesh([theta1,0],[theta2,0],nt,ratiot);
theta=theta(:,1)';




a=zeros(nr*nt,2);
for i=1:nr
  
  d=norm(apts(i,:)-bpts(i,:));
  alpha=atan2(bpts(i,2)-apts(i,2),apts(i,1)-bpts(i,1));
  beta=pi-2*alpha;
  
  r=d*sin(alpha)/sin(beta)
  c=
  
  for j=1:nt
    a((j-1)*nr+i,:)=r*[cos(theta(j)),sin(theta(j))];
  end
  
end
