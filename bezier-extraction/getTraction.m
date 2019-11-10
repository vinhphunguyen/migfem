function [traction, tangent] = getTraction (jump,gp)

global damage damage0 strength dummy  Gc

delta0   = strength / dummy;
deltaF   = 2.0 * Gc / strength;
dFmd0Inv = 1. / ( deltaF - delta0 );

tangent    = zeros(2,2);
traction   = zeros(2,1);

tension = 1;

jum     = jump(2);

if jum < 0
    tension = 0;
end

delta   = jum * tension;
hisDam  = ( delta - delta0 ) * dFmd0Inv;
damMax  = damage0(gp);

if hisDam > damMax
    loading = 1;
else
    loading = 0;
    hisDam  = damMax;
end

if hisDam > 0
    dam     = hisDam * deltaF / delta;
else
    dam    = 0.;
    hisDam = 0;
end

secant = (1-dam*tension)*dummy;

tangent(2,2)  = secant;     % be careful, with this
tangent(1,1)  = dummy;
traction(2)   = secant*jum;
traction(1)   = dummy*jump(1);

if ( loading )        
   H  = deltaF * delta0 * dFmd0Inv / ( delta * delta * delta );    
   tangent(2,2) = tangent(2,2) - dummy * H * jum * jum;
end

damage(gp) = hisDam;

