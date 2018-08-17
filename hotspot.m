function [rhoad_s_o, tsstoo]=hotspot(hot, ala, K_em, k_em, psi, w_em, rhoad_s_o, ts, to)

% Marie Weiss
% INRA Avignon

% Taking into account the ellipsoidal distribution of  mean inclination angles
% of leaves in the hot spot term.
% The HotSpot is implemented according to the changes made by
% Andrieu, B., Baret, F., Jacquemoud, S., Malthus, T. and Steven, M., 1997.
% Evaluation of an improved inversion of SAIL model for simulating bidirectional reflectance of sugar beet canopies. Remote Sensing of Environment 60, 247-257.

ko= K_em;
ks= k_em;
tgs=tan(ts);
tgo=tan(to);
tss=exp(-ks);
rsod=rhoad_s_o;
w=w_em;

dso=sqrt(abs(tgs.^2+tgo.^2-2*tgs.*tgo.*cos(psi)));

% Computation of hot-spot term;
sl = hot*pi/4/(1+0.357*(ala/(97-ala))^1.252).*sqrt(ko.*ks); % Andrieu 1997 (Eq.4)
alf= 1e6;
if hot>0
    alf=dso./sl;
end
sumint=zeros(size(alf));

tsstoo = zeros(size(tss));
ind=find(alf==0);

tsstoo(ind) = tss(ind);
sumint(ind)= (1-tss(ind))./ks(ind);
clear ind

ind=find(alf~=0);
if ~isempty(ind)
    fhot(ind)=(ko(ind).*ks(ind)).^0.5;
    x1=zeros(size(to(ind)));
    y1=zeros(size(to(ind)));
    f1=ones(size(to(ind)));
    fint=(1-exp(-alf)).*.05;

    for istep=1:20
        if istep < 20
            x2=-log(1-istep.*fint(ind))./alf(ind);
        else
            x2=ones(size(ind));
        end
        y2=-(ko(ind)+ks(ind)).*x2+fhot(ind)'.*(1-exp(-alf(ind).*x2))./alf(ind); %'
        f2=exp(y2);
        sumint(ind)=sumint(ind)+(f2-f1).*(x2-x1)./(y2-y1);
        x1=x2;
        y1=y2;
        f1=f2;
    end
    tsstoo(ind)=f1;
end

rsos=w.*sumint;
rso=rsos+rsod;

rhoad_s_o=rso;
