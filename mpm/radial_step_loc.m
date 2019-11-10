function [sigma,C,alpha,epsilonp] = radial_step_loc(epsilon,epsilonp0,alpha0,mu,kappa,k1,k0)
% performs a radial return step at one Gauss point

eye2   = [1 1 0]';
eye2x2 = [1 1 0; 
          1 1 0; 
          0 0 0];
I_dev  = eye(3) - 0.5*eye2x2;
I      = [1 0 0; 0 1 0; 0 0 0.5];
Iinv   = [1 0 0; 0 1 0; 0 0 2];
C_el   = kappa*eye2x2 + 2*mu*I*I_dev; 

% Compute trial stress
epsilon_dev  = I_dev*epsilon;
s_trial      = 2*mu*I*(epsilon_dev-epsilonp0);
norm_s_trial = sqrt(s_trial(1)^2 + s_trial(2)^2 + 2*s_trial(3)^2);
sigma_trial  = kappa*sum(epsilon(1:2))*eye2 + s_trial;

% Check yield condition
f_trial = norm_s_trial - (k1*alpha0 + k0);

if f_trial <= 0 % elastic step       
    alpha    = alpha0;
    epsilonp = epsilonp0;
    C        = C_el;
    sigma    = sigma_trial;
else % plastic step
    normal = s_trial/norm_s_trial;
    lambda = (norm_s_trial - k1*alpha0 - k0)/(2*mu + k1);
    alpha = alpha0 + lambda;
    % Update plastic strain and stress
    epsilonp = epsilonp0 + lambda*Iinv*normal; 
    sigma = kappa*sum(epsilon(1:2))*eye2 + s_trial - 2*mu*lambda*normal;
    theta1 = 1 - 2*mu*lambda/norm_s_trial;
    theta2 = 2*mu/(2*mu+k1) - 1 + theta1;
    C = kappa*eye2x2 + 2*mu*theta1*I*I_dev - 2*mu*theta2*normal*normal';
end
