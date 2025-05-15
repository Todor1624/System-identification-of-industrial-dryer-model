clearvars
close all
clc






%% Zadatak 1А: Snimanje odziva sistema na step promenu upravljanja oko nominalne radne tačke.

% Definisanje parametara eksperimenta
u_step_vent = 3;
t_step_vent = 200;
t_step_nom = 15;
u_nom_vent = 6;
u_nom_temp = 4.5;
sim_step = 0.01;
%%

% Snimljeni su rezultati eksperimenta perturbovanja sistema oko nominalne
% radne tačke.
load("izlazna_merenja.mat")


% Izdvajanje rezultata simulacije
t_sim = out.y_out_vent.Time;
t_end = t_sim(end);

y_vent = out.y_out_vent.Data;
u_vent = out.u_out_vent.Data;


% Snimljeni podaci su veoma zašumljeni. Potrebno je filtrirati ih pre dalje
% primene.

y_filt_vent = doFilter(y_vent, 0.985);

figure(1)

plot(t_sim,y_vent, 'Color','r'); hold all;
plot(t_sim,y_filt_vent, 'Color','blue');
plot(t_sim,u_vent, 'Color','black', LineStyle='--')

xlim([180 220]);
ylim([4 9.5]);
grid on

xlabel("Vreme [s]");
ylabel("Protok vazduha u tunelu [V]")
legend("Nefiltrirano merenje", "Filtrirano merenje", "Upravljanje",...
    'Location','northwest')
    
title(["Prikaz zadate sekvence upravljanja"; "i odziva sistema (STEP pobuda)"])


%% Zadatak 1B: Određivanje First Order plus Dead Time (FOPDT )aproksimacije sistema. 

% Dinamika sistema estimirana je na osnovu usponskog vremena, pojačanja u
% stacionarnom stanju.


% Pojačanje se estimira kao odnos srednje vrednosti signala u stacionarnom
% stanju i intenziteta pobude.

y_0 = mean(y_filt_vent(t_sim < 200));
y_end = mean(y_filt_vent(t_sim>210));

u_0 = u_nom_vent;
u_end = u_nom_vent + u_step_vent;

delta_y = y_end - y_0;
delta_u = u_end - u_0;

Kp = delta_y/delta_u; 

% Vizualnom inspekciom ne primećuje se značajno transportno kašnjenje, tako
% da se ono neće uračunati u aproksimaciju.

% Usponsko vreme se računa kao vreme od zadavanja pobude (pretpostavlja se da nema transportnog kašnjenja)
% do dostizanja 63% vrednosti signala u stacionarnom stanju.

y_63 = y_0 + delta_y*0.63;

t_0 = find(t_sim > 200, 1, 'first');
t_63 = find(y_filt_vent(t_0:end) > y_63, 1, 'first');
t_63 = t_0 + t_63;

% U sistemima sa jednim polom se dominantna vremenska konstanta može uzeti
% ista kao vreme uspona.

t_63 = t_sim(t_63) - t_sim(t_0); 

s = tf('s');
G_fopdt = Kp/(s*t_63 + 1);

% Perioda odabiranja se bira tako da obuhvati 20 odbiraka signala tokom
% usponskog vremena (konzervativno).


T_samp = t_63/20;     


figure(2)

u_model = u_step_vent * ones(length(y_filt_vent(t_0:end)), 1);
y_model = lsim(G_fopdt,u_model,t_sim(t_0:end));

plot(t_sim(t_0:end), y_model + y_0, "Color","blue", linewidth = 1.2); hold all
plot(t_sim(t_0:end), y_filt_vent(t_0:end), "Color","red")
grid on;
xlim([200 t_end]);

xlabel("Vreme [s]");
ylabel("Protok vazduha u tunelu [V]")
legend("Odziv FOPDT modela", "Filtriran odziv merenja",...
    'Location','southeast')
title(" Poređenje odziva filtriranog merenja i FOPDT modela")

%% Zadatak 1C: Pododabiranje merenja procenjenom periodom odabiranja.

% Odlučeno je izvršiti pododabiranje merenja, budući da je izabrana perioda
% odabiranja i više nego dovoljna da obuhvati ključnu dinamiku sistema.

clc


T_samp = t_63/20;                  

t_resampled = 0:T_samp:t_sim(end);

N_resampled = ceil(length(y_filt_vent)/length(t_resampled));

y_resampled = y_filt_vent(1:N_resampled:end);
t_resampled = t_sim(1:N_resampled:end); 
u_resampled= u_vent(1:N_resampled:end);
T_samp = t_resampled(2) - t_resampled(1);

% Odnos dužina vektora merenja i vektora uniformnih vremenskih odbiraka
% mora biti, tako da se radi blaga korekcija periode odabiranja. Odnos
% dominantne vremenske konstante i perioda odabiranja je promenjen na
% 17.85, što je i dalje više nego dovoljno.

figure
plot(t_sim, y_filt_vent,'Color', 'blue'); hold all
plot(t_resampled, y_resampled, 'Color', 'magenta', Marker='*');

xlim([150, 152]);
ylim([5.3, 5.7]);

xlabel("Vreme [s]");
ylabel("Protok vazduha u tunelu [V]")
legend("Snimljeno u 25001 tačaka", "Pododabiranje sa 3572 tačaka",...
    'Location','northeast')
title("Zumirani prikaz filtriranih merenja pre i nakon pododabiranja")
grid on;

%% Zadatak 1D: Nalaženje ARX1 modela diskretizacijom FOPDT aproksimacije sistema.



a1 = exp(pole(G_fopdt)*T_samp);
z = tf('z', T_samp);

G_fopdt_arx1 = Kp*(1-a1)/(z - a1);


y_model = lsim(u_resampled,G_fopdt,t_resampled); 
y_model_arx1 = lsim(u_resampled,G_fopdt_arx1,t_resampled); 

figure
subplot(2,1,1)
plot(t_resampled,y_model,LineWidth=1.2,LineStyle="--", Color='r'); hold all
stairs(t_resampled,y_model_arx1,LineWidth=1, Color='blue');
xlim([200 210]);
ylabel("Protok vazduha u tunelu [V]")
title("Uticaj diskretizacije na aproksimaciju prvog reda")
grid on;
legend("Kontinualni FOPDT model", "ARX1 model", 'Location','southeast')

subplot(2,1,2)
plot(t_resampled,y_model,LineWidth=1.2,LineStyle="--", Color='r'); hold all
stairs(t_resampled,y_model_arx1,LineWidth=1, Color='blue');
xlim([200 201]);
ylabel("Protok vazduha u tunelu [V]")
xlabel("Vreme [s]")
xlabel("")
grid on;

%% Zadatak 2А: Snimanje odziva na pobudne signale velikog reda perzistencije.

%% Pseudonasumična binarna sekvenca (PRBS)

T_prekidanje = 5;
N_reg = 10;
u_max = 3;
t_end = 400;
t_prbs = 250;

load("PRBS_merenje.mat");

% U ovom zadatku, cilj je identifikacija sistema koji je opisan protokom u tunelu
% industrijske sušare koristeći ARX modele prvog i drugog reda. Prvi korak je akvizicija
% merenja protoka kada je sistem doveden u nominalno stanje i pobuđen
% signalom velikog stepena perzistencije. Prvi signal koji će se koristiti
% za ovu namenu je PRBS. 

% Parametri ovakvog pobudnog signala su perioda prekidanja, trajanje, broj
% prekidačkih registara i amplituda. Prva dva parametara određuju dužinu
% sekvence i zajedno sa trećim diktiraju koji je frekvencijski opseg
% validne estimacije. Amplituda se bira tako da se sistem održi u uskom
% opsegu oko nominalne radne tačke.

% Trajanje je zadato a trajanje periode se bira tako da je obuhvaćeno vreme
% uspona i smirenja odziva. Kombinacijama registara koji učestvuju u zbiru
% po modulu 2 se postižu različiti nivoi perzistencije. Za potrebe ovog
% zadatka PRBS ostvaruje dovoljan nivo perzistencije u većini slučajeva,
% budući da se ne estimira više od tri parametara. 

% Opis sprovedenog eksperimenta je sledeći:
% 1. Dovođenje sistema u nominalan radan režim [0 - 200] s.
% 2. Zadavanje step perturbacije inteziteta 3 V [200 - 250] s.
% 3. Zadavanje unipolarne PRBS perturbacije amplitude 3V [250 - 400] s.

% Snimljeni podaci se prikazuju filtrirani, na isti način kao merenja u
% zadatku 1.

figure
plot(t_prbs_vent,u_prbs_vent, Color='black', LineStyle="--"); hold all
plot(t_prbs_vent,y_prbs_vent, Color='r', LineWidth=1.2);
grid on;
ylabel("Protok vazduha u tunelu [V]")
ylim([4.5 9.5])
xlabel("Vreme [s]")
title(["Prikaz zadate sekvence upravljanja"; "i odziva sistema (PRBS pobuda)"])
legend(["Upravljačka sekvenca", "Odziv sistema"], 'Location', 'northwest')

figure
plot(t_prbs_temp,u_prbs_temp, Color='black', LineStyle="--"); hold all
plot(t_prbs_temp,y_prbs_temp, Color='r', LineWidth=1.2);
grid on;
ylabel("Temperatura vazduha u tunelu [V]")
ylim([2.5 5])
xlabel("Vreme [s]")
title("Odziv temperature na zadato nominalno upravljanje")
legend(["Nominalno upravljanje", "Odziv sistema"], 'Location', 'northwest')


%% Bipolarna povorka (BP)

T_prekidanje = 10;
u_max = 3;
t_end = 400;
t_bp = 250;

load('BP_merenje.mat')

% Cilj ovog eksperimenta se podudara sa prethodnim. U ovom slučaju se
% umesto PRBS pobude koristi pobuda četvrtkama koje su istih širina.
% Bira se manje parametara u poređenju sa PRBS-om (perioda prekidanja i
% amplituda) a pristup pri izboru je identičan.

% Eksperiment se izvodi isto kao sa PRBS pobudom sa filtriranjem merenja.

figure
plot(t_bp_vent,u_bp_vent, Color='black', LineStyle="--"); hold all
plot(t_bp_vent,y_bp_vent, Color='r', LineWidth=1.2);
grid on;
ylabel("Protok vazduha u tunelu [V]")
ylim([4.5 9.5])
xlabel("Vreme [s]")
title(["Prikaz zadate sekvence upravljanja"; "i odziva sistema (BIPOLARNA pobuda)"])
legend(["Upravljačka sekvenca", "Odziv sistema"], 'Location', 'northwest')

figure
plot(t_bp_temp,u_bp_temp, Color='black', LineStyle="--"); hold all
plot(t_bp_temp,y_bp_temp, Color='r', LineWidth=1.2);
grid on;
ylabel("Temperatura vazduha u tunelu [V]")
ylim([2.5 5])
xlabel("Vreme [s]")
title("Odziv temperature na zadato nominalno upravljanje")
legend(["Nominalno upravljanje", "Odziv sistema"], 'Location', 'northwest')


%% Zadatak 2B: Priprema snimljenih podataka za proces identifikacije

% U narednim koracima se vrši dalja obrada snimljenih podataka, nakon
% izvršenog filtriranja. 

% Izdvajanje dela odziva na pobudu visokog reda perzistencije (samo PRBS i BP)

load("PRBS_merenje.mat");
load("BP_merenje.mat");

t_180 = find(t_prbs_vent>180,1,'first');
t_prbs_vent = t_prbs_vent(t_180:end);
y_prbs_vent = y_prbs_vent(t_180:end);
u_prbs_vent = u_prbs_vent(t_180:end);


t_180 = find(t_bp_vent>180,1,'first');
t_bp_vent = t_bp_vent(t_180:end);
y_bp_vent = y_bp_vent(t_180:end);
u_bp_vent = u_bp_vent(t_180:end);

% Detrending (izbacivanje srednjih srednjih vrednosti)

y_bp_vent = y_bp_vent - mean(y_bp_vent);
u_bp_vent = u_bp_vent - mean(u_bp_vent);

y_prbs_vent = y_prbs_vent - mean(y_prbs_vent);
u_prbs_vent = u_prbs_vent - mean(u_prbs_vent);

% Procena kašnjenja korišćenjem kroskorelacije za izolovan deo merenja

t_220 = find(t_prbs_vent>220,1,'first');
tau = length(1:t_220)/2;
rxy_prbs = kroskorelacija(y_prbs_vent(1:t_220),u_prbs_vent(1:t_220),tau);
rxy_bp = kroskorelacija(y_bp_vent(1:t_220),u_bp_vent(1:t_220),tau);
tau_arr = -tau:tau;

tau_max_prbs = tau_arr(find(rxy_prbs==max(rxy_prbs),1,'first'));
tau_max_bp = tau_arr(find(rxy_bp==max(rxy_bp),1,'first'));

figure
subplot(2,1,1)
plot(t_prbs_vent(1:t_220),u_prbs_vent(1:t_220),LineStyle="--",Color="black"); hold all
plot(t_prbs_vent(1:t_220),y_prbs_vent(1:t_220),Color="blue");
grid on
xlim([180 220])
xlabel("Vreme [s]");
ylabel("Protok vazduha u tunelu [V]")
title(["Ocena transportnog kašnjenja na osnovu kroskorelacije";"(pre PRBS pobude)"])

subplot(2,1,2)
plot(tau_arr, rxy_prbs, Color="blue");
grid on
xlabel("Kašnjenje [odb]");
ylabel("Kroskorelacija ulaz / izlaz");
xlim([-100 100])
ylim([0.4 1])

text(40,1.1,strcat("Kašnjenje [s] = ", num2str(tau_max_prbs*T_samp)))



figure
subplot(2,1,1)
plot(t_bp_vent(1:t_220),u_bp_vent(1:t_220),LineStyle="--",Color="black"); hold all
plot(t_bp_vent(1:t_220),y_bp_vent(1:t_220),Color="blue");
grid on
xlim([180 220])
xlabel("Vreme [s]");
ylabel("Protok vazduha u tunelu [V]")
title(["Ocena transportnog kašnjenja na osnovu kroskorelacije";"(pre BP pobude)"])

subplot(2,1,2)
plot(tau_arr, rxy_bp, Color="blue");
grid on
xlabel("Kašnjenje [odb]");
ylabel("Kroskorelacija ulaz / izlaz");
xlim([-100 100])
ylim([0.4 0.9])

text(40,0.98,strcat("Kašnjenje [s] = ", num2str(tau_max_bp*T_samp)))

%%

% Procenjena kašnjenja se uklanjaju iz identifikacionih objekata.

y_prbs_id = y_prbs_vent(tau_max_prbs:end);
u_prbs_id = u_prbs_vent(1:end-tau_max_prbs+1);
y_bp_id = y_bp_vent(tau_max_bp:end);
u_bp_id = u_bp_vent(1:end-tau_max_bp+1);

t_prbs_id = t_prbs_vent(1:end-tau_max_prbs+1);
t_bp_id = t_bp_vent(1:end-tau_max_bp+1);

% Na kraju se od obrađenih merenja odziva i upravljanja pravi ID objekat
id_data_prbs = iddata(y_prbs_id,u_prbs_id,T_samp);
id_data_bp = iddata(y_bp_id,u_bp_id,T_samp);

figure
subplot(2,1,1)
plot(t_prbs_id,u_prbs_id,LineStyle="--",Color="black"); hold all
plot(t_prbs_id,y_prbs_id,Color="blue");
xlim([t_prbs_id(1) t_prbs_id(end)])
ylim([-2 3])
grid on
ylabel("Protok vazduha u tunelu [V]")
legend(["Upravljačka sekvenca (PRBS)", "Odziv sistema"], 'Location', 'northwest')
title("Konačne vremenske serije podataka za identifikaciju")
subplot(2,1,2)

plot(t_bp_id,u_bp_id,LineStyle="--",Color="black"); hold all
plot(t_bp_id,y_bp_id,Color="blue");
xlim([t_prbs_id(1) t_prbs_id(end)])
ylim([-2 3])
grid on
ylabel("Protok vazduha u tunelu [V]")
legend(["Upravljačka sekvenca (BP)", "Odziv sistema"], 'Location', 'northwest')

%% Zadatak 2C: 

% id_data_prbs == podaci za identifikaciju sa PRBS signalom
% id_data_bp == podaci za identifikaciju sa BP signalom

% t_prbs_id == vremenski niz za id_data_prbs
% t_bp_id == vremenski niz za id_data_bp

% *System Identification Toolbox*    

% Koriscenjem funkcija system identification toolbox-a 
% odredili smo ARX modele prvog i drugog reda i export-ovali ih u workspace
%% PRBS identifikacija
systemIdentification('ARX_PRBS');
% !!! Exportovati modele u workspace(prevlacenje modela na polje ToWorkspace)

%% Nule i polovi DT ARX modela prvog i drugog reda

nuleARX_221_prbs=zero(arx221_prbs);
poloviARX_221_prbs=pole(arx221_prbs);

nuleARX_111_prbs=zero(arx111_prbs);
poloviARX_111_prbs=pole(arx111_prbs);
%% Nule i polovi CT ARX modela prvog i drugog reda

ARX111_CT_prbs=d2c(arx111_prbs);  % konverzija iz DT u CT domen
ARX221_CT_prbs=d2c(arx221_prbs);

% Odredjivanje nula i polova

poloviARX_111_CT_prbs = pole(ARX111_CT_prbs);
nuleARX_111_CT_prbs= zero(ARX111_CT_prbs);

poloviARX_221_CT_prbs = pole(ARX221_CT_prbs);
nuleARX_221_CT_prbs = zero(ARX221_CT_prbs);
% Iscrtavanje step odziva DT i CT modela na osnovu kojih odredjujemo
% vremenske konstante i pojacanje modela (moze i stepinfo)
figure
step(arx111_prbs);   % DT ARX model prvog reda
figure
step(ARX111_CT_prbs); % CT ARX model prvog reda

figure
step(arx221_prbs);   % DT ARX model drugog reda
figure
step(ARX221_CT_prbs); % CT ARX model drugog reda


%% BP identifikacija
systemIdentification('ARX_BP');
% !!! Exportovati modele u workspace(prevlacenje modela na polje ToWorkspace)
%% Nule i polovi DT ARX modela prvog i drugog reda
% Odredjivanje nula i polova DT modela
nuleARX_221_bp=zero(arx221_bp);
poloviARX_221_bp=pole(arx221_bp);
nuleARX_111_bp=zero(arx111_bp);
poloviARX_111_bp=pole(arx111_bp);
%% Nule i polovi CT ARX modela prvog i drugog reda
ARX111_CT_bp=d2c(arx111_bp);  % konverzija iz DT u CT domen
ARX221_CT_bp=d2c(arx221_bp);

% Odredjivanje nula i polova
poloviARX_111_CT_bp = pole(ARX111_CT_bp);
nuleARX_111_CT_bp= zero(ARX111_CT_bp);
poloviARX_221_CT_bp = pole(ARX221_CT_bp);
nuleARX_221_CT_bp = zero(ARX221_CT_bp);

% Iscrtavanje step odziva DT i CT modela na osnovu kojih odredjujemo
% vremenske konstante i pojacanje modela (moze i stepinfo)
figure
step(arx111_bp);   % DT ARX model prvog reda
figure
step(ARX111_CT_bp); % CT ARX model prvog reda

figure
step(arx221_bp);   % DT ARX model drugog reda
figure
step(ARX221_CT_bp); % CT ARX model drugog reda


%% Identifikacija WRLLS metodom

%% WRLLS PRBS identifikacija

%% ARX model prvog reda za ulazni PRBS signal 
N=length(t_prbs_id);
for i=1:3
MatrixRo=[0.95 2 10]; % faktor zaboravljanja
ro=MatrixRo(i);
Teta(:,1) = [0; 0];
P = eye(2)*100; % P=c^2*I c>>1
for k = 2:N
    fi(:,k) = [-y_prbs_id(k-1) u_prbs_id(k-1)]; %trenutni regresioni vektor za ARX
    e=y_prbs_id(k)-fi(:,k)'*Teta(:,k-1);
    P=1/ro*(P-P*fi(:,k)*fi(:,k)'*P/(ro+fi(:,k)'*P*fi(:,k)));
    K=P*fi(:,k);
    Teta(:,k)=Teta(:,k-1)+K*e;
end


a1 = Teta(1,:);
b1 =  Teta(2,:);



p_est=-a1;   % diskretni pol
K_est=b1./(1+a1);  % procena pojacanja
sp_est=-log(p_est)/T_samp;

G_est1=tf(K_est(end)*sp_est(end), [1 sp_est(end)]); % CT fja prenosa
figure
subplot(2,1,1)
step(G_est1);     % Step odziv na osnovu koga odredjujemo pojacanje i vremenske konstante
subplot(2,1,2)
pzmap(G_est1);    % polovi i nule CT funkcije prenosa

y_est =  lsim(G_est1,u_prbs_id,t_prbs_id);
figure
plot(t_prbs_id,y_prbs_id);
hold on
plot(t_prbs_id,y_est,'r--');
hold off
xlim([179 401]);
title(['Prikaz odziva stvarnog sistema i ARX modela prvog reda za ro= ' num2str(ro)]);
legend('original', 'ARX1');
end

% Na osnovu dobijenih grafika moze se uociti da izborom faktora
% zaboravljanja koji je blizak 1, algoritam ce losije pratiti brze promene
% parametara, dok za vrednosti koje su vise udaljene od 1 algoritam ce
% bolje pratiti brze promene parametara, ali losije ce pratiti sporo
% promenljive parametre.

%% ARX model drugog reda za ulazni PRBS signal

N=length(t_prbs_id);
for i=1:3
ro=MatrixRo(i);   % faktor zaboravljanja
Teta2(:,1) = zeros(4,1);
P = eye(4)*100; % P=c^2*I c>>1
for k = 3:N
    fi2(:,k) = [-y_prbs_id(k-1) -y_prbs_id(k-2) u_prbs_id(k-1) u_prbs_id(k-2)]; %trenutni regresioni vektor za ARX
    e=y_prbs_id(k)-fi2(:,k)'*Teta2(:,k-2);
    P=1/ro*(P-P*fi2(:,k)*fi2(:,k)'*P/(ro+fi2(:,k)'*P*fi2(:,k)));
    K=P*fi2(:,k);
    Teta2(:,k-1)=Teta2(:,k-2)+K*e;
end



a1=Teta2(1,end);
a2=Teta2(2,end);
b1=Teta2(3,end);
b2=Teta2(4,end);
G_ztf = tf([b1 b2],[1 a1 a2],T_samp);
G_est2=d2c(G_ztf);


figure
subplot(2,1,1)
step(G_est2);     % Step odziv na osnovu koga odredjujemo pojacanje i vremenske konstante
subplot(2,1,2)
pzmap(G_est2);    % polovi i nule CT funkcije prenosa

y_est =  lsim(G_est2,u_prbs_id,t_prbs_id);
figure
plot(t_prbs_id,y_prbs_id);
hold on
plot(t_prbs_id,y_est,'g--');
hold off
xlim([179 401]);
title(['Prikaz odziva stvarnog sistema i ARX modela drugog reda za ro= ' num2str(ro)]);
legend('original', 'ARX2');
end

% Na osnovu dobijenih grafika moze se uociti da izborom faktora
% zaboravljanja koji je blizak 1, algoritam ce losije pratiti brze promene
% parametara, dok za vrednosti koje su vise udaljene od 1 algoritam ce
% bolje pratiti brze promene parametara, ali losije ce pratiti sporo
% promenljive parametre.

%% Poredjenje odziva stvarnog sistema, ARX modela prvog i modela drugog reda na jednom grafiku za ro=10

% Izabrali smo ro koje odgovara dinamici naseg sistema, radi boljeg
% pracenja brzih promena parametara

y_est1 =  lsim(G_est1,u_prbs_id,t_prbs_id);   % estimirani odzivi sistema modeliranih kao ARX modeli prvog i drugog reda
y_est2 =  lsim(G_est2,u_prbs_id,t_prbs_id);

figure
plot(t_prbs_id,y_prbs_id);
hold on
plot(t_prbs_id,y_est1,'r--');
plot(t_prbs_id,y_est2,'g--');
hold off
xlim([179 401]);
title('Poredjenje odziva na PRBS stvarnog sistema, ARX1 i ARX2 modela')
legend('original','ARX1','ARX2')

% Sagledavanje prelaznih rezima
figure
plot(t_prbs_id,y_prbs_id);
hold on
plot(t_prbs_id,y_est1,'r--');
plot(t_prbs_id,y_est2,'g--');
hold off
xlim([250 300]);
title('Sagledavanje prelaznih rezima PRBS')
legend('original','ARX1','ARX2')


%% WRLLS BP identifikacija

%% ARX model prvog reda za ulaznu povorku bipolarnih cetvrtki 
N=length(t_bp_id);
for i=1:3
ro=MatrixRo(i); % faktor zaboravljanja
Teta11(:,1) = [0; 0];
P = eye(2)*100; % P=c^2*I c>>1
for k = 2:N
    fi11(:,k) = [-y_bp_id(k-1) u_bp_id(k-1)]; %trenutni regresioni vektor za ARX
    e=y_bp_id(k)-fi11(:,k)'*Teta11(:,k-1);
    P=1/ro*(P-P*fi11(:,k)*fi11(:,k)'*P/(ro+fi11(:,k)'*P*fi11(:,k)));
    K=P*fi11(:,k);
    Teta11(:,k)=Teta11(:,k-1)+K*e;
end


a1 = Teta11(1,:);
b1 =  Teta11(2,:);



p_est=-a1;   % diskretni pol
K_est=b1./(1+a1);  % procena pojacanja
sp_est=-log(p_est)/T_samp;

G_est11=tf(K_est(end)*sp_est(end), [1 sp_est(end)]); % CT fja prenosa
figure
subplot(2,1,1)
step(G_est11);     % Step odziv na osnovu koga odredjujemo pojacanje i vremenske konstante
subplot(2,1,2)
pzmap(G_est11);    % polovi i nule CT funkcije prenosa

y_est =  lsim(G_est11,u_bp_id,t_bp_id);
figure
plot(t_bp_id,y_bp_id);
hold on
plot(t_bp_id,y_est,'r--');
hold off
xlim([179 401]);
title(['Prikaz odziva stvarnog sistema i ARX modela prvog reda za ro= ' num2str(ro)]);
legend('original', 'ARX1');
end

% Na osnovu dobijenih grafika moze se uociti da izborom faktora
% zaboravljanja koji je blizak 1, algoritam ce losije pratiti brze promene
% parametara, dok za vrednosti koje su vise udaljene od 1 algoritam ce
% bolje pratiti brze promene parametara, ali losije ce pratiti sporo
% promenljive parametre.

%% ARX model drugog reda za ulaznu povorku bipolarnih cetvrtki

N=length(t_bp_id);
for i=1:3
ro=MatrixRo(i);   % faktor zaboravljanja
Teta22(:,1) = zeros(4,1);
P = eye(4)*100; % P=c^2*I c>>1
for k = 3:N
    fi22(:,k) = [-y_bp_id(k-1) -y_bp_id(k-2) u_bp_id(k-1) u_bp_id(k-2)]; %trenutni regresioni vektor za ARX
    e=y_bp_id(k)-fi22(:,k)'*Teta22(:,k-2);
    P=1/ro*(P-P*fi22(:,k)*fi22(:,k)'*P/(ro+fi22(:,k)'*P*fi22(:,k)));
    K=P*fi22(:,k);
    Teta22(:,k-1)=Teta22(:,k-2)+K*e;
end



a1=Teta22(1,end);
a2=Teta22(2,end);
b1=Teta22(3,end);
b2=Teta22(4,end);
G_ztf = tf([b1 b2],[1 a1 a2],T_samp);
G_est22=d2c(G_ztf);


figure
subplot(2,1,1)
step(G_est22);     % Step odziv na osnovu koga odredjujemo pojacanje i vremenske konstante
subplot(2,1,2)
pzmap(G_est22);    % polovi i nule CT funkcije prenosa

y_est =  lsim(G_est22,u_bp_id,t_bp_id);
figure
plot(t_bp_id,y_bp_id);
hold on
plot(t_bp_id,y_est,'g--');
hold off
xlim([179 401]);
title(['Prikaz odziva stvarnog sistema i ARX modela drugog reda za ro= ' num2str(ro)]);
legend('original', 'ARX2');
end

% Na osnovu dobijenih grafika moze se uociti da izborom faktora
% zaboravljanja koji je blizak 1, algoritam ce losije pratiti brze promene
% parametara, dok za vrednosti koje su vise udaljene od 1 algoritam ce
% bolje pratiti brze promene parametara, ali losije ce pratiti sporo
% promenljive parametre.

%% Poredjenje odziva stvarnog sistema, ARX modela prvog i modela drugog reda na jednom grafiku za ro=10

% Izabrali smo ro koje odgovara dinamici naseg sistema, radi boljeg
% pracenja brzih promena parametara

y_est1 =  lsim(G_est11,u_bp_id,t_bp_id);   % estimirani odzivi sistema modeliranih kao ARX modeli prvog i drugog reda
y_est2 =  lsim(G_est22,u_bp_id,t_bp_id);

figure
plot(t_bp_id,y_bp_id);
hold on
plot(t_bp_id,y_est1,'r--');
plot(t_bp_id,y_est2,'g--');
hold off
xlim([179 401]);
title('Poredjenje odziva na BP stvarnog sistema, ARX1 i ARX2 modela')
legend('original','ARX1','ARX2')

% Sagledavanje prelaznih rezima
figure
plot(t_bp_id,y_bp_id);
hold on
plot(t_bp_id,y_est1,'r--');
plot(t_bp_id,y_est2,'g--');
hold off
xlim([250 300]);
title('Sagledavanje prelaznih rezima BP')
legend('original','ARX1','ARX2')






