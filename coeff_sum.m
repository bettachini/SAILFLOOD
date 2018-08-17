function  [k_sum,s_sum,a_sum,sigma_sum,s_prima_sum,w_sum,v_sum,u_sum,K_sum]=coeff_sum(h,ts_w,tl,aw,bw,L_prime_sum,refl,tran,to_w,psi, theta_crit)



% María Eugenia Beget1 and Víctor Bettachini2 
% 1 Los Reseros y Las Cabañas s/n, 1686 Hurlingham, Buenos Aires, Argentina. Tel.: +54 11 4621 1684; fax: +54 11 4621 5663.
% mbeget@cnia.inta.gov.ar 
% 2 Laboratorio de Optoelectrónica, Instituto Tecnológico de Buenos Aires, Av. Eduardo Madero 399, C1106ACD Buenos Aires, Argentina.
% vbettach@itba.edu.ar

% Beget, M.E., V. A. Bettachini, C. Di Bella, F. Baret. 2013.  
% SAILHFlood: a radiative transfer model for flooded vegetation. Ecological
% Modelling, 257:25-35.

% Computation of SAIL coefficients for submerged vegetation layer, modified from
% Verhoef, W., 1984. Light scattering by leaf layers with application to canopy reflectance modeling: the SAIL model. Remote Sensing of Environment 16, 125-141.

format long e;

% k_sum: canopy extintion + water (scattering + absorption) 
if (ts_w+tl)<(pi/2)                                 % Verhoef 1984 p.131
    beta_s_w= pi;
else
    beta_s_w= acos(-1/(tan(ts_w)*tan(tl)));         % Verhoef 1984 (Eq.9)
end
k_leaf_sum= (2* L_prime_sum/ pi)* ((beta_s_w-(pi/2))*cos(tl)+ sin(beta_s_w)*tan(ts_w)*sin(tl));          % SAIL coeff Verhoef 1984 (Eq.22)
k_water= (1/ cos(ts_w))* (aw+bw);	% water (scattering + absorption) integrated for water path	% Beget et al. 2013 (Eq.21)
k_sum= k_leaf_sum +  k_water;	% SAIL coeff + water (scattering + absorption)	% Beget et al. 2013 (Eq.22)

% Angular distribution of Rayleigh scattering calculated for a water
% molecule is axially symmetric and the probabilities of forward and backward
% scattering are equal.  (http://www.wetlabs.com/iopdescript/scatterimages/fig12.gif)

s_leaf_sum= ((refl + tran)/2) *  k_leaf_sum - L_prime_sum * ((refl - tran)/2) * cos(tl)^2;              % SAIL Verhoef 1984 (Eq.30)
s_water= (1/cos(ts_w ) )* bw/ 2;	% solar flux downward, direct -> diffuse, 1/2 of scattering coefficients for pure water % Beget et al. 2013, eq. 23
s_sum= s_leaf_sum+ s_water;	% SAIL + water scattering	% Beget et al. 2013 (Eq.24)

s_prime_leaf_sum= ((refl + tran)/2) *  k_leaf_sum + L_prime_sum * ((refl - tran)/2) * cos(tl)^2;        % SAIL Verhoef 1984 (Eq.29)
s_prime_water= (1/ cos(ts_w ) )* bw/ 2;              % solar flux upward, direct -> diffuse, 1/2 of scattering coefficients for pure water % Beget et al. 2013, eq. 23
s_prima_sum= s_prime_leaf_sum+ s_prime_water;   % SAIL + water scattering	% Beget et al. 2013 (Eq.24)

% Forward scattering 
sigma_leaf_sum = L_prime_sum * ( ((refl + tran)/2) + ((refl - tran)/2) * cos(tl)^2 );                   % SAIL Verhoef 1984 (Eq.26)
if (bw+aw)==0
    sigma_sum= sigma_leaf_sum;
else
    int_sigma_water= @(thetaw_s) (h/ (2*(1- theta_crit) ) )* bw* tan(thetaw_s );
    sigma_water= quad(int_sigma_water, 0, theta_crit);	% Beget et. al 2013 (Eq.25)

%     iteraciones=20;                             % sigma_water is calculated with an average of 20 angles for the computation of path because is for diffuse flux (all directions)  
%     sigma_water=0;
%     jump=h*(aw+bw);
%     theta_crit=acos(jump);
%     if (theta_crit>pi/2)
%         theta_crit=pi/2;
%     end
%     for i=0:iteraciones-1
%         sigma_water= sigma_water + (1/iteraciones) * (h * (bw/2)/( cos((i/iteraciones)* theta_crit)));
%     end
    sigma_sum= sigma_leaf_sum + sigma_water;     % SAIL + water scattering	% Beget et al. 2013 (Eq.26)
end


% Backward scattering 
sigma_prime_leaf_sum = L_prime_sum * ( ((refl + tran)/2) - ((refl - tran)/2) * cos(tl)^2 );             % SAIL Verhoef 1984 (Eq.27)
if (bw+aw)==0
    sigma_prima_sum= sigma_prime_leaf_sum;
else
   sigma_prime_water= sigma_water;	% Beget et. al 2013 (Eq.25) 

%     iteraciones=20;
%     sigma_prime_water=0;
%     jump=h*(aw+bw);
%     theta_crit=acos(jump);
%     if (theta_crit>pi/2)
%         theta_crit=pi/2;
%     end
%     for i=0:iteraciones-1
%         sigma_prime_water= sigma_prime_water + (1/iteraciones) * (h * (bw/2)/( cos((i/iteraciones)* theta_crit) ));
%     end
    sigma_prima_sum = sigma_prime_leaf_sum + sigma_prime_water;  % SAIL + water scattering
end


a_leaf_sum = L_prime_sum - sigma_prime_leaf_sum;                % SAIL Verhoef 1984 (Eq.24) 
if (bw+aw)==0
    a_sum= a_leaf_sum;
else
    int_a_water= @(thetaw_s) (h/(1- cos(theta_crit) ))* (aw+ 0.5* bw)* tan(thetaw_s);
    a_water= quad( int_a_water, 0, theta_crit);	% Beget et. al 2013 (Eq.27) 

%     iteraciones=20;
%     a_water=0;
%     jump=h*(aw+bw);
%     theta_crit=acos(jump);
%     if (theta_crit>pi/2)
%         theta_crit=pi/2;
%     end
%     for i=0:iteraciones-1
%         a_water= a_water + (1/iteraciones) * (jump /( cos((i/iteraciones)*theta_crit) ));
%     end
%     a_water= a_water - sigma_prime_water;                          % to not to compute the diffuse forward scattering twice, because a fraction continues its way in the same direction via forward scattering, and thus does not contribute to attenuation                         
   a_sum= a_leaf_sum + a_water;	% SAIL Verhoef84 (24) +  water forward scattering and absorption	% Beget et. al 2013, eq. 28
end

% beta_o
if (to_w+tl)<(pi/2)                             % Verhoef 1984 p.131
    beta_o_w= pi;
else
    beta_o_w= acos(-1/(tan(to_w)*tan(tl)));     % Verhoef 1984 (Eq.9)
end

K_leaf_sum=(2/pi)*L_prime_sum*((beta_o_w - (pi/2)) * cos(tl) + sin(beta_o_w) * tan(to_w) * sin(tl));               % SAIL Verhoef 1984 (Eq.23)
K_water= (1/ cos(to_w ) ) * (aw+ bw );	% Beget et al. 2013, eq. 29
K_sum=K_leaf_sum +  K_water;	% SAIL coeff + water (scattering + absorption)	% Beget et al. 2013 (Eq.31)

v_sum= ((refl + tran)/2) * K_sum + ((refl - tran)/2) * L_prime_sum * cos(tl)^2;	% Verhoef 1984 (Eq.31)	% Beget et al. 2013 (Eq.32)

% proposal:
% v_hoja_sum= ((refl + tran)/2) * K_leaf_sum + ((refl - tran)/2) * L_prime_sum * cos(tl)^2; 
% if (bw+aw)==0
%     v_sum= v_hoja_sum;
% else
%     iteraciones=20;                             %v_agua is calculated with an average of 20 angles for the computation of path because is for diffuse flux (all directions)  
%     v_agua=0;
%     for i=0:iteraciones-1
%         v_agua= v_agua + (1/iteraciones) * (h * (bw/2)/( cos((i/iteraciones)*(pi/2))));
%     end
%     v_sum= v_hoja_sum + v_agua;     % SAIL + water scattering
% end   

u_sum= ((refl + tran)/2) * K_sum - ((refl - tran)/2) * L_prime_sum * cos(tl)^2;	% Verhoef84 (32)	% Beget et al. 2013 (Eq.32)

% proposal:
% u_hoja_sum= ((refl + tran)/2) * K_sum - ((refl - tran)/2) * L_prime_sum *  cos(tl)^2;                % Verhoef 1984 (Eq.32)
% if (bw+aw)==0
%     u_sum= u_hoja_sum;
% else
%     iteraciones=20;                             %u_agua is calculated with an average of 20 angles for the computation of path because is for diffuse flux (all directions)  
%     u_agua=0;
%     for i=0:iteraciones-1
%         u_agua= u_agua + (1/iteraciones) * (h * (bw/2)/( cos((i/iteraciones)*(pi/2))));
%     end
%     u_sum= u_hoja_sum + u_agua;     % SAIL + water scattering
% end   

% Auxiliary azimuthal angles, Verhoef 1984 p.133
if (psi<=abs(beta_s_w - beta_o_w))
    beta1_w= psi;
    beta2_w= abs(beta_s_w - beta_o_w);
    beta3_w= 2*pi - beta_s_w - beta_o_w;
elseif (psi > abs(beta_s_w - beta_o_w) && psi < (2*pi - beta_s_w - beta_o_w))
    beta1_w= abs(beta_s_w - beta_o_w);
    beta2_w= psi;
    beta3_w= 2*pi - beta_s_w - beta_o_w;
elseif (psi>= (2*pi - beta_s_w - beta_o_w))
    beta1_w= abs(beta_s_w - beta_o_w);
    beta2_w= 2*pi - beta_s_w - beta_o_w;
    beta3_w= psi;
end

w_sum= (L_prime_sum / (2*pi)) * ((pi * refl - beta2_w * (refl + tran)) * (2*cos(tl)^2 + sin(tl)^2 * tan(ts_w) * tan(to_w) * cos(psi)) + (refl + tran) * sin(beta2_w) * (2*cos(tl)^2/(cos(beta_s_w) * cos(beta_o_w)) + cos(beta1_w) * cos(beta3_w) * sin(tl)^2 * tan(ts_w) * tan(to_w) ));	% Verhoef84 (33)	% Beget et al. 2013, eq. 33

% proposal:
% w_hoja_sum= (L_prime_sum / (2*pi)) * ((pi * refl - beta2_w * (refl + tran)) * (2*cos(tl)^2 + sin(tl)^2 * tan(ts_w) * tan(to_w) * cos(psi)) + (refl + tran) * sin(beta2_w) * (2*cos(tl)^2/(cos(beta_s_w) * cos(beta_o_w)) + cos(beta1_w) * cos(beta3_w) * sin(tl)^2 * tan(ts_w) * tan(to_w) ));                 % Verhoef84 (33)
% w_agua=(h/cos(ts_w))* (bw/2);  
% w_sum=w_hoja_sum + w_agua;
