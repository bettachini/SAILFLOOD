function  [k_em,s_em,a_em,sigma_em,s_prima_em,w_em,v_em,u_em,K_em]=coeff_em(ts,tl,L_prime_em,refl,tran,to,psi)


% María Eugenia Beget1 and Víctor Bettachini2 
% 1 Los Reseros y Las Cabañas s/n, 1686 Hurlingham, Buenos Aires, Argentina. Tel.: +54 11 4621 1684; fax: +54 11 4621 5663.
% mbeget@cnia.inta.gov.ar 
% 2 Laboratorio de Optoelectrónica, Instituto Tecnológico de Buenos Aires, Av. Eduardo Madero 399, C1106ACD Buenos Aires, Argentina.
% vbettach@itba.edu.ar

% Beget, M.E., V. A. Bettachini, C. Di Bella, F. Baret. 2013.  
% SAILHFlood: a radiative transfer model for flooded vegetation. Ecological
% Modelling, 257:25-35.

% Computation of SAIL coefficients for emerged vegetation layer from
% Verhoef, W., 1984. Light scattering by leaf layers with application to canopy reflectance modeling: the SAIL model. Remote Sensing of Environment 16, 125-141.

if (ts+tl)<(pi/2)                                                                           % Verhoef 1984 p.131
    beta_s= pi;
else
    beta_s= acos(-1/(tan(ts)*tan(tl)));                                                     % Verhoef 1984 (Eq.9)
end
k_em= (2/pi)* L_prime_em* ((beta_s - (pi/2))* cos(tl)+ sin(beta_s)* tan(ts)* sin(tl));      % atenuation Verhoef 1984 (Eq.22)

s_em= ((refl + tran)/2) *  k_em - L_prime_em * ((refl - tran)/2) * cos(tl)^2;               % Verhoef 1984 (Eq.30)
s_prima_em= ((refl + tran)/2) *  k_em + L_prime_em * ((refl - tran)/2) * cos(tl)^2;         % Verhoef 1984 (Eq.29)
sigma_em = L_prime_em * ( ((refl + tran)/2) + ((refl - tran)/2) * cos(tl)^2 );              % Verhoef 1984 (Eq.26)
sigma_prime_em = L_prime_em * ( ((refl + tran)/2) - ((refl - tran)/2) * cos(tl)^2 );        % Verhoef 1984 (Eq.27)
a_em = L_prime_em - sigma_prime_em;        

if (to+tl)<(pi/2)                                                                           % Verhoef 1984 p.131
    beta_o= pi;
else
    beta_o= acos(-1/(tan(to)*tan(tl)));                                                     % Verhoef 1984 (Eq.9)
end
K_em=(2/pi)*L_prime_em*((beta_o - (pi/2)) * cos(tl) + sin(beta_o) * tan(to) * sin(tl));     % Verhoef 1984 (Eq.23)

v_em= ((refl + tran)/2) * K_em + ((refl - tran)/2) * L_prime_em * cos(tl)^2;                % Verhoef 1984 (Eq.31)
u_em= ((refl + tran)/2) * K_em - ((refl - tran)/2) * L_prime_em * cos(tl)^2;                % Verhoef 1984 (Eq.32)

% Azimuth auxiliar angles, Verhoef 1984 p.133
if (psi<=abs(beta_s - beta_o))
    beta1= psi;
    beta2= abs(beta_s - beta_o);
    beta3= 2*pi - beta_s - beta_o;
elseif (psi > abs(beta_s - beta_o) && psi < (2*pi - beta_s - beta_o))
    beta1= abs(beta_s - beta_o);
    beta2= psi;
    beta3= 2*pi - beta_s - beta_o;
elseif (psi>= (2*pi - beta_s - beta_o))
    beta1= abs(beta_s - beta_o);
    beta2= 2*pi - beta_s - beta_o;
    beta3= psi;
end

w_em= (L_prime_em / (2*pi)) * ((pi * refl - beta2 * (refl + tran)) * (2*cos(tl)^2 + sin(tl)^2 * tan(ts) * tan(to) * cos(psi)) + (refl + tran) * sin(beta2) * (2*cos(tl)^2/(cos(beta_s) * cos(beta_o)) + cos(beta1) * cos(beta3) * sin(tl)^2 * tan(ts) * tan(to) ));                 % Verhoef84 (33)
