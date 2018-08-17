function Reflectance=run_sailhflood


% María Eugenia Beget1 and Víctor Bettachini2 
% 1 Los Reseros y Las Cabañas s/n, 1686 Hurlingham, Buenos Aires, Argentina. Tel.: +54 11 4621 1684; fax: +54 11 4621 5663.
% mbeget@cnia.inta.gov.ar 
% 2 Laboratorio de Optoelectrónica, Instituto Tecnológico de Buenos Aires, Av. Eduardo Madero 399, C1106ACD Buenos Aires, Argentina.
% vbettach@itba.edu.ar

% Beget, M.E., V. A. Bettachini, C. Di Bella, F. Baret. 2013.  
% SAILHFlood: a radiative transfer model for flooded vegetation. Ecological
% Modelling, 257:25-35.

% Run SAILHFLood model
% Subrutines required: coeff_sum.m, coeff_em.m, coeff_ab.m, coeff_nw.m,  hotspot.m, 
% phittospurom_tobira.m (Leaf optical properties), soil.m (Soil optical properties),  
% --> Leaf optical properties and soil optical properties files
% can be replaced by another plant and soil respectively. Please find and
% replace rutine name in sailhflood.m!!!

% INPUTS

% Incident radiation
to=15*(pi/180);     % INPUT
ts=60*(pi/180);     % INPUT
psi= 0* pi/180;     % INPUT
skyl= 0;            % INPUT

%Create zeros vector Reflectance
Reflectance(:,1)=[400:2400]; % lambda
Reflectance(:,2)=zeros(size(Reflectance(:,1))); 

%Five cases varying LAI and water height (h, cm)

%Case 1
LAI_sum=	0.5	;     % INPUT
LAI_em=	0.5	;     % INPUT
h=	5	;      % INPUT
for lambda=400:2400
    [Reflectance(lambda-399,2)]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda);
end
R1=Reflectance;

%Case 2
LAI_sum=	1	;     % INPUT
LAI_em=	1	;     % INPUT
h=	1	;     % INPUT
for lambda=400:2400
    [Reflectance(lambda-399,2)]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda);
end
R2=Reflectance;

%Case 3
LAI_sum=	1.5	;     % INPUT
LAI_em=	1.5	;     % INPUT
h=	2	;     % INPUT
for lambda=400:2400
    [Reflectance(lambda-399,2)]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda);
end
R3=Reflectance;

%Case 4
LAI_sum=	2	;     % INPUT
LAI_em=	2	;     % INPUT
h=	4	;     % INPUT

for lambda=400:2400
    [Reflectance(lambda-399,2)]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda);
end
R4=Reflectance;

%Case 5
LAI_sum=	2.5	;     % INPUT
LAI_em=	2.5	;     % INPUT
h=	5	;     % INPUT
for lambda=400:2400
    [Reflectance(lambda-399,2)]=sailh_flood(LAI_sum,LAI_em,h,ts,to,psi,skyl,lambda);
end
R5=Reflectance;

% save run_sailhflood % save workspace as run_sailhflood.m
%---------------------------------------------------
figure('Color',[1 1 1]);hold all;
plot(Reflectance(:,1),R1,'-r');
plot(Reflectance(:,1),R2,'-g');
plot(Reflectance(:,1),R3,'-b');
plot(Reflectance(:,1),R4,'-m');
plot(Reflectance(:,1),R5,'-c');
axis([400 2400 0 1]);
xlabel('Wavelength (nm)');set(gca,'XTick',400:200:2400);
ylabel('Reflectance');set(gca,'YTick',0:0.2:1);
box('on');legend('R1','R2','R3','R4','R5');
            
            
            
