function [Refl_total]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda)

% María Eugenia Beget1 and Víctor Bettachini2 
% 1 Los Reseros y Las Cabañas s/n, 1686 Hurlingham, Buenos Aires, Argentina. Tel.: +54 11 4621 1684; fax: +54 11 4621 5663.
% mbeget@cnia.inta.gov.ar 
% 2 Laboratorio de Optoelectrónica, Instituto Tecnológico de Buenos Aires, Av. Eduardo Madero 399, C1106ACD Buenos Aires, Argentina.
% vbettach@itba.edu.ar

% Beget, M.E., V. A. Bettachini, C. Di Bella, F. Baret. 2013.  
% SAILHFlood: a radiative transfer model for flooded vegetation. Ecological
% Modelling, 257:25-35.

% sailh_flood.m:  SAILH model modified to compute reflectance of a canopy partcially submerged by
% pure water. % Some parts of code are taked from sail2.m (Marie Weiss)
% SAIL model
% Verhoef, W., 1984. Light scattering by leaf layers with application to canopy reflectance modeling: the SAIL model. Remote Sensing of Environment 16, 125-141.
% Verhoef, W., 1985. Earth observation modeling based on layer scattering matrices. Remote Sensing of Environment 17, 165-178.
% Andrieu, B., Baret, F., Jacquemoud, S., Malthus, T. and Steven, M., 1997. Evaluation of an improved inversion of SAIL model for simulating bidirectional reflectance of sugar beet canopies. Remote Sensing of Environment 60, 247-257.


format long e;
% INPUTS
% LAI_sum: submerged  LAI
% LAI_em: emerged LAI
% h: water height (m)
% ts: solar zenith angle (degrees)
% to: view zenith angle (degrees)
% psi: azimutal angle difference between solar and view angles (degrees)
% skyl: diffuse fraction
% lambda: wavelength (nm)

% Subrutines required: coeff_sum.m, coeff_em.m, coeff_ab.m, coeff_nw.m,  hotspot.m, 
% phittospurom_tobira.m (Leaf optical properties), soil.m (Soil optical properties),  
% --> Leaf optical properties and soil optical properties files
% can be replaced by another plant and soil respectively. Please find and replace rutine name in this code!!!

% SAIL INPUTS that are fixed, but can be MODIFIED as inputs in other
% subrutines (by replacing leaf_optical_properties.m and
% soil_optical_properties.m)
% ala: mean leaf inclination angle (phittospurom_tobira.m)
% hot: hot spot parameter (phittospurom_tobira.m)
% refl: leaf reflectance (phittospurom_tobira.m)  
% tran: leaf transmittance (phittospurom_tobira.m)
% Soil_r_s_o, Soil_r_s_d, Soil_r_d_o, Soil_r_d_d: soil reflectance
% coefficients  (soil.m)

% Water optical properties (absoption and scattering coefficients)are
% computed in coeff_ab.m
% Water refractive coefficient is computed in coeff_nw

% OUTPUT
% [Refl_total]: total bidirectional reflectance 


% Angles: degrees to radian
to= to *(pi/180);  % to: observation zenith angle (0,pi)
ts= ts *(pi/180);  % ts: sola zenith angle (0,pi)
psi= psi *(pi/180); % observation azimuth angle - solar azimuth angle (0,2*pi)
if to<0
    psi= psi+pi;
    if psi> 2*pi
        psi= psi-2*pi;
    end
end

% Azimuthal angle
if to<0
    psi=psi+pi;
    if psi>2*pi
        psi=psi-2*pi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leaf inclination distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ala, hot, refl, tran]= phittospurom_tobira(lambda );

excent=exp(-1.6184e-5.*ala.^3+2.1145e-3.*ala.^2-1.2390e-1.*ala+3.2491);
tlref=[5 15 25 35 45 55 64 72 77 79 81 83 85 87 89];
alpha1=[10 20 30 40 50 60 68 76 78 80 82 84 86 88 90];
alpha2=[0 alpha1(1,1:14)];
ff=zeros(15,0);

% -----------relative frequencies calculation

tlref=tlref'*pi/180;          
alpha1=alpha1'*pi/180;         
alpha2=alpha2'*pi/180;         
x1=excent./sqrt(1+excent.^2..*tan(alpha1).^2);
x2=excent./sqrt(1+excent.^2..*tan(alpha2).^2);
if excent == 1
    ff=abs(cos(alpha1)-cos(alpha2));
else
    a1=excent./sqrt(abs(1-excent.^2));
    a12=a1.^2;
    x12=x1.^2;
    x22=x2.^2;
    a1px1=sqrt(a12+x12);
    a1px2=sqrt(a12+x22);
    a1mx1=sqrt(a12-x12);
    a1mx2=sqrt(a12-x22);
    if excent >1
        ff=x1.*a1px1+a12.*log(x1+a1px1);
        ff=abs(ff-(x2.*a1px2+a12.*log(x2+a1px2)));
    else
        ff=x1.*a1mx1+a12.*asin(x1./a1);
        ff=abs(ff-(x2.*a1mx2+a12.*asin(x2./a1)));
    end
end
ff=ff./sum(ff);


% NOTATION
% Ed_c_m
% d= direction of propagation (u: upward,d: downward)
% m= medium (w: water, a: air)
% c= vertical dimension (0: top of canopy, 1: water-air interface , 2: soil)
%             d         u
% |Es|   | tss tds | rds ros |  |Es|
% |E-|   | tsd tdd | rdd rod |  |E-|
%     = --------------------- 
% |E+|   | rsd rdd | tdd tod |  |E+|
% |Eo|   | rso rdo | tdo too |  |Eo|
%             d         u
% REFLECTANCE MATRICES
% R_c_m
% c= vertical dimension (0: top of canopy, 1: water-air interface , 2: soil)
% m= medium (w: water, a: air)
% TRANSMITTANCE MATRICES
% T_c_m
% c= vertical dimension (0: top of canopy, 1: water-air interface , 2: soil)
% m= medium (w: water, a: air)
% REFLECTANCE MATRICES ELEMENTS
% rhowd_i_s
% d= direction of propagation (u: upward,d: downward)
% i= INCIDENT (s: solar, d: diffuse)
% s= OUTWARD (d: difusse, o: observer)
% TRANSMITTANCE MATRICES ELEMENTS
% taud_i_s
% d= direction of propagation (u: upward,d: downward)
% i= INCIDENT (s: solar, d: diffuse)
% s= OUTWARD (d: difusse, o: observer)
% SOIL REFLECTANCE MATRIX ELEMENTS
% r_i_s
% i= INCIDENT (s: solar, d: diffuse)
% s= OUTWARD (d: difusse, o: observer)
% COEFFICIENTS
% coef_em: SAIL coeff for emerged vegetation
% coef_sum: final coeff for submerged vegetation (coef_leaf_sum + coef_water)
% coef_leaf_sum:  SAIL coeff for submerged vegetation 
% coef_water: pure water coeff

if h==0             % WITHOUT WATER

    % DRY SOIL
    % SOIL REFLECTANCE MATRIX
    [r_s_d, r_d_d, r_s_o, r_d_o]= soil(lambda, ts, to);
    R_2_s(1,1)=r_s_d;   % directional-hemisferical reflectance, input
    R_2_s(1,2)=r_d_d;   % bihemisferical reflectance, input
    R_2_s(2,1)=r_s_o;   % bidirectional reflectance, input (s->o)
    R_2_s(2,2)=r_d_o;   % hemi-directional reflectance, 
    
    R_11_iws=R_2_s;     % if there is no water -> the soil matrix is assigned as interface reflectance matrix
else

% SOIL + WATER

n_w=coeff_nw(lambda);        % WATER REFRACTIVE INDEX
[aw,bw]=coeff_ab(lambda);    % aw (lambda) Absorption coefficient (aw) m-1     
                           % bw (lambda) Scattering coefficient (bw) m-1

% Critical angle (total internal reflection)
theta_crit= asin((1/ n_w ) );	% Beget et. al 2013, (Eq.17)


% Change in ts due to water interface (ts_w) (Snell)
ts_w= asin(sin(ts)/n_w);     % zenital solar angle below water interface

% Change in to due to water interface (to_w) (Snell)
to_w= asin( sin(to ) /n_w );     % zenital observation angle below water interface	% Beget et al. 2013 (Eq.30)

% WET SOIL
% SOIL REFLECTANCE MATRIX
[r_s_d, r_d_d, r_s_o, r_d_o]=soil(lambda,ts,to);  
R_2_s(1,1)=r_s_d;   % directional-hemisferical reflectance, input
R_2_s(1,2)=r_d_d;   % bihemisferical reflectance, input
R_2_s(2,1)=r_s_o;   % bidirectional reflectance, input (s->o)
R_2_s(2,2)=r_d_o;   % hemi-directional reflectance, 

% SUBMERGED VEGETATION COEFFICIENTS
k_sum=0;
s_sum=0;
a_sum=0;
sigma_sum=0;
s_prime_sum=0;
w_sum=0;
v_sum=0;
u_sum=0;
K_sum=0;
for i=1:length(tlref)
    tl=tlref(i);
    L_prime_sum=LAI_sum*ff(i);
    [k_sumP,s_sumP,a_sumP,sigma_sumP,s_prime_sumP,w_sumP,v_sumP,u_sumP,K_sumP]=coeff_sum(h,ts_w,tl,aw,bw,L_prime_sum,refl,tran,to_w,psi, theta_crit);  
    k_sum=k_sum+k_sumP;
    s_sum=s_sum+s_sumP;
    a_sum=a_sum+a_sumP;
    sigma_sum=sigma_sum+sigma_sumP;
    s_prime_sum=s_prime_sum+s_prime_sumP;
    w_sum=w_sum+w_sumP;
    v_sum=v_sum+v_sumP;
    u_sum=u_sum+u_sumP;
    K_sum=K_sum+K_sumP;
end

m_sum = sqrt(a_sum^2 - sigma_sum^2);                % Verhoef 1985 p.167
h1_sum = (a_sum + m_sum) / (sigma_sum);             % Verhoef 1985 (Eq.9a)
h2_sum = (a_sum - m_sum) / (sigma_sum);             % Verhoef 1985 (Eq.9b)

% SUBMERGED MATRICES

%             d         u
% |Es|   | tss tds | rds ros |  |Es|
% |E-|   | tsd tdd | rdd rod |  |E-|
%     = --------------------- 
% |E+|   | rsd rdd | tdd tod |  |E+|
% |Eo|   | rso rdo | tdo too |  |Eo|
%             d         u

% Downward AND upward FLUXES ARE ASSUMED SYMMETRICAL

% Downward
% leaves + water(scattering + absorption)
tauwd_s_s=exp(- ( k_sum ));   % Verhoef 1985 (Eq.15) + k_sum

if exp(m_sum)>1.771849450036987e+305    % a_sum takes a high value => exp(m_sum) = Inf; it uses a value previos to crash
    tauwd_d_d=0;
    rhowd_d_d=1/h1_sum;
else
    tauwd_d_d= (h1_sum - h2_sum)/(h1_sum * exp(m_sum) - h2_sum * exp(-m_sum));                      % Verhoef 1985 (Eq.15)
    rhowd_d_d=(exp(m_sum)-exp(-m_sum))/(h1_sum * exp(m_sum) - h2_sum * exp(-m_sum));                % Verhoef 1985 (Eq.15)
end
%as sigma_sum=0: (there is no water nor leafs to disperse)
% h1_sum = INF
% h2_sum = NAN

D_s_sum= (-s_sum * (k_sum + a_sum) - s_prime_sum * sigma_sum) /( (k_sum)^2 - (m_sum)^2);        % Verhoef 1985 (Eq.11a)
C_s_sum= (s_prime_sum * ((k_sum - a_sum)) - s_sum * sigma_sum)/( (k_sum)^2 - (m_sum)^2);        % Verhoef 1985 (Eq.10a)

tauwd_s_d= D_s_sum * (tauwd_s_s - tauwd_d_d) - C_s_sum * tauwd_s_s * rhowd_d_d;                 % Verhoef 1985 (Eq.15)
rhowd_s_d= C_s_sum * (1 - tauwd_s_s * tauwd_d_d) - D_s_sum * rhowd_d_d;                         % Verhoef 1985 (Eq.15)
tauwu_o_o= exp(- ( K_sum ));                                                                    % Verhoef 1985 (Eq.15)
tauwu_d_d= tauwd_d_d;                                                                           % ARE ASSUMED SYMMETRICAL

D_o_sum= (-u_sum*(K_sum+a_sum)- v_sum * sigma_sum)/(K_sum^2 - m_sum^2);                         % Verhoef 1985 (Eq.11b)
C_o_sum= (v_sum*(K_sum - a_sum) - u_sum * sigma_sum)/(K_sum^2 - m_sum^2);                       % Verhoef 1985 (Eq.10b)

rhowd_d_o= C_o_sum * (1 - tauwu_o_o * tauwd_d_d) - D_o_sum * rhowd_d_d;                         % Verhoef 1985 (Eq.15) 

if hot>0 % HOTSPOT
    H_o_sum= (s_sum * C_o_sum + s_prime_sum * D_o_sum)/(k_sum + K_sum);                         % Verhoef 1985 HOTSPOT
else
    H_o_sum= (s_sum * C_o_sum + s_prime_sum * D_o_sum + w_sum)/(k_sum + K_sum);                 % Verhoef 1985 (Eq.12b)
end
rhowd_s_o= H_o_sum * (1 - tauwd_s_s * tauwu_o_o) - C_o_sum * tauwd_s_d * tauwu_o_o - D_o_sum * rhowd_s_d;   % Verhoef 1985 (Eq.15)
rhowu_d_d= rhowd_d_d;                                                                           % ARE ASSUMED SYMMETRICAL               
tauwu_d_o= D_o_sum * (tauwu_o_o - tauwu_d_d) - C_o_sum * tauwu_o_o * rhowu_d_d;                 % Verhoef 1985 (Eq.15)
tauwd_d_s= 0;                                                                                   % Verhoef 1984 (Eq.1a) : Suits does not consider another source of "s" than the sun: we assume intensity of diffuse radiation is lower than solar radiation
tauwu_o_d= 0;                                                                                   % Verhoef 1984 (Eq.1b) : Suits does not consider diffuse radiation in observer direction: doble scattering 

% Transmitance matrix for water layer (downward fluxes)
T_12_w(1,1)=tauwd_s_s;
T_12_w(1,2)=tauwd_d_s;          % = 0 
T_12_w(2,1)=tauwd_s_d;
T_12_w(2,2)=tauwd_d_d;

% Transmitance matrix for water layer (upward fluxes)
T_2_w(1,1)=tauwu_d_d;
T_2_w(1,2)=tauwu_d_o;           % = 0 
T_2_w(2,1)=tauwu_o_d;
T_2_w(2,2)=tauwu_o_o;

% Reflectance matrix for water layer (downward fluxes, top)
R_12_w(1,1)=rhowd_s_d;
R_12_w(1,2)=rhowd_d_d;
R_12_w(2,1)=rhowd_s_o;
R_12_w(2,2)=rhowd_d_o;

% Reflectance matrix for water layer (upward fluxes, bottom)
                        
R_2_w(1,1)=0;           % ds = 0 , doble scattering is not contemplated
R_2_w(1,2)=0;           % os = 0 , doble scattering is not contemplated
R_2_w(2,1)=rhowu_d_d;
R_2_w(2,2)=0;           % od = 0 , doble scattering is not contemplated

% REFLECTANCE MATRIX, SOIL + WATER

R_12_ws= T_2_w * R_2_s * (eye(2) - R_2_w * R_2_s )^(-1) * T_12_w +  R_12_w; % Verhoef 1985:  R^*_t


% HOTSPOT
if hot>0
    [rhowd_s_o, tsstoo]= hotspot(hot, ala, K_sum, k_sum, psi, w_sum, rhowd_s_o, ts, to);
    R_12_ws(2,1)= rhowd_s_o+ tsstoo* r_s_o+ (1/ (1- r_d_d* rhowu_d_d ) )* ( ...
    (tauwd_s_s * r_s_d + tauwd_s_d * r_d_d) * tauwu_d_o + ...
    (tauwd_s_d + tauwd_s_s * r_s_d * rhowu_d_d) * r_d_o * tauwu_o_o );
end

% INTERFACE + (SOIL + WATER)
%             d         u
% |Es|   | tss tds | rds ros |  |Es|
% |E-|   | tsd tdd | rdd rod |  |E-|
%     = --------------------- 
% |E+|   | rsd rdd | tdd tod |  |E+|
% |Eo|   | rso rdo | tdo too |  |Eo|
%             d         u

% DOWNWARD FLUXES INTERFACE COEFFICIENTS
rho_descendente= @(theta).5* ((sin(theta- asin(sin(theta )/ n_w) )./sin(theta+ asin(sin(theta )/ n_w) ) ).^2+ (tan(theta- asin(sin(theta )/ n_w ) )./tan(theta+ asin(sin(theta )/ n_w) ) ).^2 );	% Beget et. al 2013, eq. 15
rhoid_d_d= quad(rho_descendente, 0, pi/2 );	% Fresnel (integration  between 0 - pi/2)

aux_refl= 5*((sin(ts-ts_w)./sin(ts+ts_w)).^2+(tan(ts-ts_w)./tan(ts+ts_w) ).^2);    % Beget et. al 2013 (Eq.13)
if ts==to
    rhoid_s_o= aux_refl;	% Beget et. al 2013 (Eq.13)
    rhoid_s_d= 0;                       % interface does not diffract

else
    rhoid_s_o= 0;
    rhoid_s_d= aux_refl; 	 % Beget et. al 2013 (Eq.13)
end

rhoid_d_o= 0;
tauid_s_s= 1- aux_refl;    % tauid_s_s: complement from reflexion: Fresnel	% Beget et. al 2013 (Eq. 14)
tauid_d_s=0;                                                                            % interface does not cause difussion
tauid_s_d=0;                                                                            % interface does not cause difussion
tauid_d_d= 1- rhoid_d_d;                                                                % complement from reflexion: Fresnel (between 0 - pi/2)	% Beget et. al 2013 (Eq.16)


 % UPWARD FLUXES INTERFACE COEFFICIENTS
rho_ascendente=@(thetaW).5* (abs(((sin(thetaW- asin(n_w* sin(thetaW) ) ) ./sin(thetaW+ asin(n_w* sin(thetaW ) ) ) ) ) ).^2+ abs((tan(thetaW- asin(n_w* sin(thetaW) ) )./tan(thetaW+ asin(n_w* sin(thetaW ) ) ) ) ).^2);
rhoiu_d_d=quad(rho_ascendente, 0, pi/2);	% Fresnel (integration  between 0 - pi/2) % Beget et al. 2013 (Eq.19)


    
tauiu_d_d= 1- rhoiu_d_d;	% Beget et al. 2013 (Eq.20)
rhoiu_d_s=0; % negligible  
rhoiu_o_s=0; % except specular
tauiu_o_d=0; % does not  scattering
tauiu_d_o=0; % complement from reflexion: Fresnel 

if to~=0
    if to_w>=theta_crit                             % total incident energy is reflected by total internal reflection (Hetch, E., 2000. Óptica. Pearson (Ed.). 720 pp.)
        rhoiu_o_d= 1;
    else
        rhoiu_o_d =(.5*((sin(to_w- to)./sin(to_w+ to ) ).^2+ (tan(to_w- to)./tan(to_w+ to) ).^2) );	% Beget et. al 2013 (Eq.18)
    end
else
    rhoiu_o_d= (n_w^2+ 1) /(1+ n_w)^2; % Si to=0 entonces to_w=0. Hecht (4.47) indice refraccion mayor a menor sumamos como difuso lo reflejado hacia abajo : no lo consideramos doble rebote
end
tauiu_o_o= 1-rhoiu_o_d;	% Beget et al. 2013, eq. 20



% Reflectance matrix for interface layer (downward fluxes, top)
R_11_i(1,1)=rhoid_s_d;
R_11_i(1,2)=rhoid_d_d;
R_11_i(2,1)=rhoid_s_o;
R_11_i(2,2)=rhoid_d_o;

% Reflectance matrix for interface layer (upward fluxes, bottom)
R_12_i(1,1)=rhoiu_d_s;      
R_12_i(1,2)=rhoiu_o_s;     
R_12_i(2,1)=rhoiu_d_d;
R_12_i(2,2)=rhoiu_o_d;     

% T_11_i % Transmitance matrix for interface layer (downward fluxes)
T_11_i(1,1)=tauid_s_s;
T_11_i(1,2)=tauid_d_s;
T_11_i(2,1)=tauid_s_d;
T_11_i(2,2)=tauid_d_d;

% T_12_i % Transmitance matrix for interface layer (upward fluxes)
T_12_i(1,1)=tauiu_d_d;
T_12_i(1,2)=tauiu_o_d;
T_12_i(2,1)=tauiu_d_o;
T_12_i(2,2)=tauiu_o_o;

% REFLECTANCE MATRIX, INTERFACE + (SOIL + WATER)

R_11_iws= T_12_i * R_12_ws * (eye(2) - R_12_i * R_12_ws )^(-1) * T_11_i +  R_11_i;

end   

if LAI_em==0
    Refl_total=real(R_11_iws(2,1)* (1- skyl )+ R_11_iws(2,2)* ( skyl ) );
else

% AIR + (INTERFACE + WATER + SOIL)
% EMERGED VEGETATION COEFFICIENTS
k_em=0;
s_em=0;
a_em=0;
sigma_em=0;
s_prime_em=0;
w_em=0;
v_em=0;
u_em=0;
K_em=0;
for i=1:length(tlref)
    tl=tlref(i);
    L_prime_em=LAI_em*ff(i);
    [k_emP,s_emP,a_emP,sigma_emP,s_prime_emP,w_emP,v_emP,u_emP,K_emP]=coeff_em(ts,tl,L_prime_em,refl,tran,to,psi);  
    k_em=k_em+k_emP;
    s_em=s_em+s_emP;
    a_em=a_em+a_emP;
    sigma_em=sigma_em+sigma_emP;
    s_prime_em=s_prime_em+s_prime_emP;
    w_em=w_em+w_emP;
    v_em=v_em+v_emP;
    u_em=u_em+u_emP;
    K_em=K_em+K_emP;
end

m_em = sqrt(a_em^2 - sigma_em^2);               % Verhoef 1985 p.167
h1_em = (a_em + m_em) / (sigma_em);             % Verhoef 1985 (Eq.9a)
h2_em = (a_em - m_em) / (sigma_em);             % Verhoef 1985 (Eq.9b)


% EMERGED VEGETATION MATRICES

% SAIL COEFF
tauad_s_s=exp(- ( k_em ));                                                                      % Verhoef 1985 (Eq.15) 
tauad_d_d= (h1_em - h2_em)/(h1_em * exp(m_em) - h2_em * exp(-m_em));                            % Verhoef 1985 (Eq.15)
rhoad_d_d=(exp(m_em)-exp(-m_em))/(h1_em * exp(m_em) - h2_em * exp(-m_em));                      % Verhoef 1985 (Eq.15)
D_s_em= (-s_em * (k_em + a_em) - s_prime_em * sigma_em) /( (k_em)^2 - (m_em)^2);                % Verhoef 1985 (Eq.11a)
C_s_em= (s_prime_em * ((k_em - a_em)) - s_em * sigma_em)/( (k_em)^2 - (m_em)^2);                % Verhoef 1985 (Eq.10a)
tauad_s_d= D_s_em * (tauad_s_s - tauad_d_d) - C_s_em * tauad_s_s * rhoad_d_d;                   % Verhoef 1985 (Eq.15)
rhoad_s_d= C_s_em * (1 - tauad_s_s * tauad_d_d) - D_s_em * rhoad_d_d;                           % Verhoef 1985 (Eq.15)
tauau_o_o= exp(- ( K_em ));                                                                     % Verhoef 1985 (Eq.15)
tauau_d_d= tauad_d_d;
D_o_em= (-u_em*(K_em+a_em)- v_em * sigma_em)/(K_em^2 - m_em^2);                                 % Verhoef 1985 (Eq.11b)
C_o_em= (v_em*(K_em - a_em) - u_em * sigma_em)/(K_em^2 - m_em^2);                               % Verhoef 1985 (Eq.10b)
rhoad_d_o= C_o_em * (1 - tauau_o_o * tauad_d_d) - D_o_em * rhoad_d_d;                           % Verhoef 1985 (Eq.15) 

if hot>0    % HOTSPOT
    H_o_em= (s_em * C_o_em + s_prime_em * D_o_em)/(k_em + K_em);                                % Verhoef 1985 (Eq.12b) HOT
else
    H_o_em= (s_em * C_o_em + s_prime_em * D_o_em + w_em)/(k_em + K_em);                         % Verhoef 1985 (Eq.12b)
tomamos el valor anterior a estallar 
end

rhoad_s_o= H_o_em * (1 - tauad_s_s * tauau_o_o) - C_o_em * tauad_s_d * tauau_o_o - D_o_em * rhoad_s_d;   % Verhoef 1985 (Eq.15)
rhoau_d_d= rhoad_d_d;                                                            
tauau_d_o= D_o_em * (tauau_o_o - tauau_d_d) - C_o_em * tauau_o_o * rhoau_d_d;                   % Verhoef 1985 (Eq.15)



% Transmitance matrix for interface layer (downward fluxes)
T_0_a(1,1)=tauad_s_s; 
T_0_a(1,2)=0;           
T_0_a(2,1)=tauad_s_d; 
T_0_a(2,2)=tauad_d_d;

% Transmitance matrix for air layer (upward fluxes)
T_11_a(1,1)=tauau_d_d; 
T_11_a(1,2)=0;          
T_11_a(2,1)=tauau_d_o; 
T_11_a(2,2)=tauau_o_o; 

% Reflectance matrix for air layer (downward fluxes, top)
R_0_a(1,1)=rhoad_s_d;
R_0_a(1,2)=rhoad_d_d;
R_0_a(2,1)=rhoad_s_o;
R_0_a(2,2)=rhoad_d_o;

% Reflectance matrix for air layer (upward fluxes, bottom)
R_11_a(1,1)=0;           
R_11_a(1,2)=0;          
R_11_a(2,1)=rhoau_d_d;
R_11_a(2,2)=0;          

% REFLECTANCE MATRIX, AIR + (INTERFACE + WATER + SOIL)

R_0_aiws= T_11_a * R_11_iws * (eye(2) - R_11_a * R_11_iws )^(-1) * T_0_a +  R_0_a; 

% HOTSPOT
if hot>0
    [rhoad_s_o,tsstoo]=hotspot(hot,ala,K_em,k_em,psi,w_em,rhoad_s_o,ts,to);
    R_0_aiws(2,1)= rhoad_s_o+ tsstoo * R_11_iws(2,1) + (1/ (1 - R_11_iws(1,2) * rhoau_d_d)) * ( ...
        (tauad_s_s * R_11_iws(1,1) + tauad_s_d * R_11_iws(1,2)) * tauau_d_o + ...
        (tauad_s_d + tauad_s_s * R_11_iws(1,1) * rhoau_d_d) * R_11_iws(2,2) * tauau_o_o );
end
Refl_total=real(R_0_aiws(2,1)*(1-skyl)+R_0_aiws(2,2)*(skyl)); % OUTPUT: TOTAL REFLECTANCE

end
