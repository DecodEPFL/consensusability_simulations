%% Simulation code for Consensusability of linear interconnected multi-agent systems (Journal) 
%  First simulation scenario: Charging a network of clustered supercapacitors 
clear all

% System parameters

% Number of supercapacitors
N = 9;
% Number of states per supercapacitor
n = 1;
% Discretization sampling period (Forward Euler)
Ts = 1e-4;

% Electrical components - Identical local parameters
R = 5.0e3;
C = 1.0e1;

a = 1-Ts/R/C;

% Graphs: Physical Graph

Bp = [1  0  0  0  0  0;
      0  1  0  0  0  0;
     -1 -1  0  0  0  0;
      0  0  1  0  0  0;
      0  0  0  1  0  0;
      0  0  0 -1  1  0;
      0  0 -1  0 -1  0;
      0  0  0  0  0  1;
      0  0  0  0  0 -1];
  
Rij = 1e-0*(40*rand(size(Bp,2),1)+10*ones(size(Bp,2),1));
  
Wp = diag(1./Rij);

Lp = Bp*Wp*Bp';
Lpm = Ts/C*Lp;


% Graphs: Communication Graph

Bc = [1  1  0  0  0  0  0  0  0  0;
     -1  0  1  0  0  0  0  0  0  0;
      0 -1  0  1  0  0  0  0  0  0;
      0  0 -1  0  1  1  0  0  0  0;
      0  0  0  0 -1  0  1  0  0  0;
      0  0  0  0  0 -1  0  1  1  0;
      0  0  0  0  0  0 -1 -1  0  0;
      0  0  0 -1  0  0  0  0  0  1;
      0  0  0  0  0  0  0  0 -1 -1];
  
Wc = diag(ones(size(Bc,2),1));

Lc = Bc*Wc*Bc';
Lcm = Ts/C*Lc;


%% Sufficient conditions for consensusability and controller design

eigLp = eig(Lpm);

lambdapmin = eigLp(2);
lambdapmax = eigLp(end);

Deltap = lambdapmax-lambdapmin;

eigLc = eig(Lcm);

lambdacmin = eigLc(2);
lambdacmax = eigLc(end);

gammac = lambdacmax/lambdacmin;

CON11 = lambdapmin>a-1
CON12 = (gammac-1)*(1-a+lambdapmin)<gammac*(2-Deltap)

CON21 = lambdapmax<1+a
CON22 = (gammac-1)*(-1-a+lambdapmax)>gammac*(Deltap-2)

if CON11 & CON12
    Kplus_l = (-1-a+lambdapmax)/(lambdacmin);
    Kplus_u = (1-a+lambdapmin)/(lambdacmax);
    
    Kint_C1 = [max([Kplus_l,0]),Kplus_u]
    k_C1 = Kint_C1*[0.2;0.8];
end
if CON21 & CON22
    Kmin_l = (-1-a+lambdapmax)/(lambdacmax);
    Kmin_u = (1-a+lambdapmin)/(lambdacmin);

    Kint_C2 = [Kmin_l,min([Kmin_u,0])]
    k_C2 = Kint_C2*[0.005;0.995];
end

%% Simulations with CT dynamics 

Tsim = 1;

k_C2 = -200;

X = zeros(N,Tsim/Ts+1);
X(:,1) = 4+2*rand(N,1);

A = -1/R/C*eye(N)-1/C*Lp+k_C2/C*Lc;
% eig(A)

time = 0:Ts:Tsim;
for ind = 1:length(time)-1
    X(:,ind+1) = expm(A*Ts)*X(:,ind);
end


int = 'Interpreter';
lat = 'Latex';
fonts = 'Fontsize';
fs = 16;
font = 'FontName';
fn = 'Times';
line = 'LineWidth';
lw = 2;


figure(1)
set(gca, font, fn, fonts, 14);
hold on
grid on

ylabel('$$\mathbf{V_i(V)}$$', font, fn, int,lat,fonts,fs)
xlabel('Time ($$\mathbf{s}$$)', font, fn, int,lat,fonts,fs)

plot(time,X',line,lw)


legend({'$$\mathbf{V_1}$$','$$\mathbf{V_2}$$','$$\mathbf{V_3}$$','$$\mathbf{V_4}$$','$$\mathbf{V_5}$$','$$\mathbf{V_6}$$'}, font, fn, int,lat,fonts,fs);


%% Trial with different k

% krange = [-logspace(log(-Kmin_l)/log(10)+0.1,-16,1000), 0, logspace(-16, log(Kplus_u)/log(10)+1, 100)];
% krange = [linspace(-50000,Kmin_l+1000,100000), linspace(Kmin_l+1000,0,10), linspace(0, 1, 100000), linspace(1, 10000, 100)];
% eigs = zeros(length(krange),1);
% ind = 1; 
krange_m = linspace(Kmin_l*1.00001,Kmin_l*0.99999,1e5);
krange_p = linspace(Kplus_u*0,Kplus_u*5,1e5);
eigs_m = zeros(length(krange_m),1);
eigs_p = zeros(length(krange_p),1);

% for k = krange
%     eigs(ind) = max(abs(eig(a*eye(N)-Lpm+k*Lcm)));
%     ind = ind+1;
% end

ind = 1;

for k = krange_m
    eigs_m(ind) = max(abs(eig(a*eye(N)-Lpm+k*Lcm)));
    ind = ind+1;
end

ind = 1;

for k = krange_p
    eigs_p(ind) = max(abs(eig(a*eye(N)-Lpm+k*Lcm)));
    ind = ind+1;
end

kreal_min = min(krange_m(eigs_m<1))
kreal_max = max(krange_p(eigs_p<1))

% figure(2)
% set(gca, font, fn, fonts, 14);
% hold on
% grid on
% ylims = [-0.1, 1.1];
% ylim(ylims);
% 
% xlabel('$$\mathbf{k}$$', font, fn, int,lat,fonts,fs)
% ylabel('Consensusability (YES=1, NO=0)', font, fn, int,lat,fonts,fs)
% 
% plot(krange, eigs<1)
% plot([Kmin_l, Kmin_l],ylims,'r--')
% plot([Kplus_u, Kplus_u],ylims,'r--')

figure(3)
subplot(1,2,1)
set(gca, font, fn, fonts, 14);
hold on
grid on
ylims = [-0.1, 1.1];
ylim(ylims);

xlabel('$$k$$', font, fn, int,lat,fonts,fs)
ylabel({'Consensusability ', '(Yes $$=~1$$, No $$=~0$$)'}, font, fn, int,lat,fonts,fs)

plot(krange_m, eigs_m<1,line,lw)
plot([Kmin_l, Kmin_l],ylims,'r--',line,lw)

% figure(2)
subplot(1,2,2)
set(gca, font, fn, fonts, 14);
hold on
grid on
ylim(ylims);

xlabel('$$k$$', font, fn, int,lat,fonts,fs)
% ylabel('Consensusability (YES=1, NO=0)', font, fn, int,lat,fonts,fs)

plot(krange_p, eigs_p<1,line,lw)
plot([Kplus_u, Kplus_u],ylims,'r--',line,lw)

%% Extended simulations for different values of Delta_p (NOT PRETTY)

% Kints_C1 = [];
% Kints_C2 = [];
% 
% scale_Lp = [1, 10, 100, 1000, 1e4];
% 
% figure(4)
% subplot(1,2,1)
% set(gca, font, fn, fonts, 14);
% hold on
% grid on
% ylims = [-0.1, 1.1];
% ylim(ylims);
% xlabel('$$k$$', font, fn, int,lat,fonts,fs)
% ylabel({'Consensusability ', '(Yes $$=~1$$, No $$=~0$$)'}, font, fn, int,lat,fonts,fs)
% subplot(1,2,2)
% set(gca, font, fn, fonts, 14);
% hold on
% grid on
% ylim(ylims);
% xlabel('$$k$$', font, fn, int,lat,fonts,fs)
% 
% colors = ['b', 'r', 'g', 'm', 'c'];
% colind=1;
% 
% for scale = scale_Lp
%     sLp = scale*Lp;
%     sLpm = Ts/C*sLp;
%     
%     seigLp = eig(sLpm);
% 
%     slambdapmin = seigLp(2);
%     slambdapmax = seigLp(end);
% 
%     sDeltap = slambdapmax-slambdapmin;
% 
%     sCON11 = slambdapmin>a-1
%     sCON12 = (gammac-1)*(1-a+slambdapmin)<gammac*(2-sDeltap)
% 
%     sCON21 = slambdapmax<1+a
%     sCON22 = (gammac-1)*(-1-a+slambdapmax)>gammac*(sDeltap-2)
% 
%     if sCON11 & sCON12
%         sKplus_l = (-1-a+slambdapmax)/(lambdacmin);
%         sKplus_u = (1-a+slambdapmin)/(lambdacmax);
% 
%         Kints_C1(end+1,:) = [max([sKplus_l,0]),sKplus_u];
%     end
%     if sCON21 & sCON22
%         sKmin_l = (-1-a+slambdapmax)/(lambdacmax);
%         sKmin_u = (1-a+slambdapmin)/(lambdacmin);
% 
%         Kints_C2(end+1,:) = [sKmin_l,min([sKmin_u,0])];
%     end
%     
%     skrange_m = linspace(sKmin_l*1.001,sKmin_l*0.999,1e5);
%     skrange_p = linspace(sKplus_u*0,sKplus_u*5,1e5);
%     seigs_m = zeros(length(skrange_m),1);
%     seigs_p = zeros(length(skrange_p),1);
% 
%     ind = 1;
% 
%     for k = skrange_m
%         seigs_m(ind) = max(abs(eig(a*eye(N)-sLpm+k*Lcm)));
%         ind = ind+1;
%     end
% 
%     ind = 1;
% 
%     for k = skrange_p
%         seigs_p(ind) = max(abs(eig(a*eye(N)-sLpm+k*Lcm)));
%         ind = ind+1;
%     end
% 
%     skreal_min = min(skrange_m(seigs_m<1));
%     skreal_max = max(skrange_p(seigs_p<1));
%     
%     subplot(1,2,1)
%     plot(skrange_m, seigs_m<1,line,lw, 'Color', colors(colind))
%     plot([sKmin_l, sKmin_l],ylims, 'Color',colors(colind), 'LineStyle', '--',line,lw)
%     subplot(1,2,2)
%     plot(skrange_p, seigs_p<1,line,lw, 'Color', colors(colind))
%     plot([sKplus_u, sKplus_u],ylims, 'Color', colors(colind),'LineStyle','--',line,lw)
%     colind = colind+1;
% end













