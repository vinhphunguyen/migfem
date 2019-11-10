function material = createOrthotropicMaterial (theta_,rank,...
    e1,e2,e3,nu12,nu23,nu31,g12,g23,g31)

pi = 3.14159265;
c  = cos( theta_ * pi / 180.0 );
s  = sin( theta_ * pi / 180.0 );
sc = s*c;
c2 = c*c;
s2 = s*s;


if ( rank == 3 )
    transformMat_     = zeros(6,6);
    transformMatInv_  = zeros(6,6);
    materialCompMat_  = zeros(6,6);
    materialStiffMat_ = zeros(6,6);
    
    transformMat_(1,1) = c2;
    transformMat_(1,2) = s2;
    transformMat_(1,4) = 2.0 * sc;
    
    transformMat_(2,1) = s2;
    transformMat_(2,2) = c2;
    transformMat_(2,4) = - 2.0 * sc;
    
    transformMat_(3,3) = 1.0;
    
    transformMat_(4,1) = - sc;
    transformMat_(4,2) = sc;
    transformMat_(4,4) = c2 - s2;
    
    transformMat_(5,5) = c;
    transformMat_(5,6) = - s;
    
    transformMat_(6,5) = s;
    transformMat_(6,6) = c;
    
    transformMatInv_(1,1) = c2;
    transformMatInv_(1,2) = s2;
    transformMatInv_(1,4) = - 2.0 * sc;
    
    transformMatInv_(2,1) = s2;
    transformMatInv_(2,2) = c2;
    transformMatInv_(2,4) = 2.0 * sc;
    
    transformMatInv_(3,3) = 1.0;
    
    transformMatInv_(4,1) = sc;
    transformMatInv_(4,2) = - sc;
    transformMatInv_(4,4) = c2 - s2;
    
    transformMatInv_(5,5) = c;
    transformMatInv_(5,6) = s;
    
    transformMatInv_(6,5) = - s;
    transformMatInv_(6,6) = c;
    
    materialCompMat_(1,1) = 1.0 / e1;
    materialCompMat_(2,2) = 1.0 / e2;
    materialCompMat_(3,3) = 1.0 / e3;
    
    materialCompMat_(1,2) = -nu12 / e1;
    materialCompMat_(2,1) = -nu12 / e1;
    materialCompMat_(1,3) = -nu31 / e1;
    materialCompMat_(3,1) = -nu31 / e1;
    materialCompMat_(2,3) = -nu23 / e2;
    materialCompMat_(3,2) = -nu23 / e2;
    
    materialCompMat_(4,4) = 1.0 / g12;
    materialCompMat_(5,5) = 1.0 / g23;
    materialCompMat_(6,6) = 1.0 / g31;
    
    materialStiffMat_ = inv(materialCompMat_);
    
elseif (rank == 2) % PLANE STRAIN ONLY
    
    transformMat_     = zeros(3,3);
    transformMatInv_  = zeros(3,3);
    materialCompMat_  = zeros(3,3);
    
    
    transformMat_(1,1) = c2;
    transformMat_(1,2) = s2;
    transformMat_(1,3) = 2.0 * sc;

    transformMat_(2,1) = s2;
    transformMat_(2,2) = c2;
    transformMat_(2,3) = - 2.0 * sc;

    transformMat_(3,1) = - sc;
    transformMat_(3,2) = sc;
    transformMat_(3,3) = c2 - s2;

    transformMatInv_(1,1) = c2;
    transformMatInv_(1,2) = s2;
    transformMatInv_(1,3) = - 2.0 * sc;

    transformMatInv_(2,1) = s2;
    transformMatInv_(2,2) = c2;
    transformMatInv_(2,3) = 2.0 * sc;

    transformMatInv_(3,1) = sc;
    transformMatInv_(3,2) = -sc;
    transformMatInv_(3,3) = c2 - s2;

    
    materialCompMat_(1,1) = 1.0 / e1;
    materialCompMat_(2,2) = 1.0 / e2;
    materialCompMat_(3,3) = 1.0 / e3;
    
    materialCompMat_(1,2) = -nu12 / e1;
    materialCompMat_(2,1) = -nu12 / e1;
    
    materialCompMat_(1,3) = -nu31 / e1;
    materialCompMat_(3,1) = -nu31 / e1;
    materialCompMat_(2,3) =  -nu23 / e2;
    materialCompMat_(3,2) =  -nu23 / e2;
    
    materialStiffMat_ = inv(materialCompMat_);
    
    materialStiffMat_(:,  3) = 0.;
    materialStiffMat_(  3,:) = 0.;
    materialStiffMat_(  3,  3)   = g12;
    
end

material.transformMat     = transformMat_;
material.transformMatInv  = transformMatInv_;
material.transformMatInvT = transformMatInv_';
material.materialStiffMat = materialStiffMat_;
material.materialCompMat  = materialCompMat_;
material.stiffMat         = transformMatInv_ * materialStiffMat_ * material.transformMatInvT;




