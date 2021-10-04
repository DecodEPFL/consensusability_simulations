%% Simulation code for Consensusability of linear interconnected multi-agent systems (Journal) 
%  First simulation scenario: Charging a network of clustered supercapacitors 
clear all

% System parameters

% Number of DGUs
N = 9;
% Number of states per DGU
n = 3;
% Discretization sampling period (Forward Euler)
Ts = 100e-6;

% Electrical components - Identical local parameters
Rt = 2.0e-1;
Ct = 2.2e-3;
Lt = 1.8e-3;
Rl = 9;

Kpr = [-2.134,-0.163,13.553];

A = [-1/Rl/Ct     , 1/Ct         , 0        ;
     (Kpr(1)-1)/Lt,(Kpr(2)-Rt)/Lt, Kpr(3)/Lt;
     -1           ,   0          , 0        ];
 
B = [0;
     0;
     1];
 
Ap = [1/Ct, 0, 0;
      0   , 0, 0;
      0   , 0, 0];
  

Ad = eye(n)+Ts*A;
Bd = Ts*B;
Apd = Ts*Ap;


% Graphs: Physical Graph

% % Bp = [1  0  0  0  0  0;
% %       0  1  0  0  0  0;
% %      -1 -1  0  0  0  0;
% %       0  0  1  0  0  0;
% %       0  0  0  1  0  0;
% %       0  0 -1 -1  1  0;
% %       0  0  0  0 -1  0;
% %       0  0  0  0  0  1;
% %       0  0  0  0  0 -1];
%   
% Bp = [Bp'; 0 0 1 0 0 0 0 -1 0]';
% Bp = [Bp'; 0 0 0 0 0 1 0 -1 0]';
% % Bp = [Bp'; 0 1 0 -1 0 0 0 0 0]';
% 
% % Bp = [Bp'; 1 -1 0 0 0 0 0 0 0]';
% % Bp = [Bp'; 1 0 0 0 -1 0 0 0 0]';
% % Bp = [Bp'; 0 0 0 0 0 0 1 0 -1]';
% % Bp = [Bp'; 0 0 0 1 -1 0 0 0 0]';
% % Bp = [Bp'; 0 0 0 1 0 -1 0 0 0]';
% % Bp = [Bp'; 0 0 0 0 1 0 -1 0 0]';
% % Bp = [Bp'; 1 0 0 0 0 0 0 0 -1]';

% BUS TOPOLOGY: 

Bp = [1  0;
     -1  1;
      0 -1;
      0  0;
      0  0;
      0  0;
      0  0;
      0  0;
      0  0];

  
Bp = [Bp';  0  1  0 -1  0  0  0  0  0]';
Bp = [Bp';  0  0  0  1 -1  0  0  0  0]';  
Bp = [Bp';  0  0  0  1  0 -1  0  0  0]';
Bp = [Bp';  0  0  0  0  0  1 -1  0  0]';  
Bp = [Bp';  0  0  0  0  0  1  0 -1  0]';
Bp = [Bp';  0  0  0  0  0  0  0  1 -1]';
% Bp = [Bp'; 1 0 0 0 -1 0 0 0 0]';
% Bp = [Bp'; 0 0 0 0 0 0 1 0 -1]';
% Bp = [Bp'; 1 0 0 0 0 0 0 0 -1]';

% Bp = [1  0  0  0  0  0;
%       0  1  0  0  0  0;
%      -1 -1  0  0  0  0;
%       0  0  1  0  0  0;
%       0  0  0  1  0  0;
%       0  0  0 -1  1  0;
%       0  0 -1  0 -1  0;
%       0  0  0  0  0  1;
%       0  0  0  0  0 -1];
%   
% Bp = [Bp'; 0 0 1 0 0 0 0 -1 0]';
% Bp = [Bp'; 0 0 0 0 0 1 0 -1 0]';
  
Rij = 4e-0*(rand(size(Bp,2),1)+ones(size(Bp,2),1));
  
Wp = diag(1./Rij);

Lp = Bp*Wp*Bp';
% Lp = Bp*Bp';
Lpm = Ts*Lp;


% Graphs: Communication Graph

% % Ring Graph:
% Bc = [1  0  0  0  0  0  0  0 -1;
%      -1  1  0  0  0  0  0  0  0;
%       0 -1  1  0  0  0  0  0  0;
%       0  0 -1  1  0  0  0  0  0;
%       0  0  0 -1  1  0  0  0  0;
%       0  0  0  0 -1  1  0  0  0;
%       0  0  0  0  0 -1  1  0  0;
%       0  0  0  0  0  0 -1  1  0;
%       0  0  0  0  0  0  0 -1  1];
  
% % Star Graph:
% Bc = [1  1  1  1  1  1  1  1;
%      -1  0  0  0  0  0  0  0;
%       0 -1  0  0  0  0  0  0;
%       0  0 -1  0  0  0  0  0;
%       0  0  0 -1  0  0  0  0;
%       0  0  0  0 -1  0  0  0;
%       0  0  0  0  0 -1  0  0;
%       0  0  0  0  0  0 -1  0;
%       0  0  0  0  0  0  0 -1];

% Complete Graph:
Bc = [];
for v = 1:N-1
    for vv = v+1:N
        Bc(:,end+1) = zeros(N,1);
        Bc(v,end) = 1;
        Bc(vv,end) = -1;
    end
end
% Bc(:,23)=[];
% Bc(:,15)=[];
% Bc(:,32)=[];
  
Wc = diag(rand(size(Bc,2),1));
Lc = Bc*Wc*Bc';

% Lc = Lp;

% Lc = 1e-0*(eye(N)-1/N*ones(N,N));
Lcm = Ts*Lc;


%% Sufficient conditions for consensusability and controller design


eigLp = eig(Lp);
Deltap = eigLp(end)-eigLp(2)

eigLc = eig(Lcm);
gammac = max(eigLc)/min(eigLc(eigLc~=min(eigLc)));

Vs = cell(N-1,1);
invVs = cell(N-1,1);
vs = cell(N-1,1);
en = zeros(1,n);
en(end) = 1;

V = zeros((N-2)*n,(N-1)*n);
v = zeros((N-2)*n,1);

for i = 2:N
    Ai = Ad-eigLp(i)*Apd;
    Bi = eigLc(i)*Bd;
%     det(Ai)
    invMri = inv(ctrb(Ai,Bi));
%     rank(Ai)
%     Ai^n
    Moi = obsv(Ai,eye(n));
    Vs{i-1} = kron(eye(n),-en*invMri)*Moi;
    invVs{i-1} = inv(Vs{i-1});
%     kron(eye(n),-en*invMri)*Moi
    vs{i-1} = -en*invMri*Ai^n;
%     -en*invMri*Ai^n
end

for i=2:N-1
    V((i-2)*n+1:(i-1)*n,(i-2)*n+1:(i-1)*n) = Vs{i-1}';
    V((i-2)*n+1:(i-1)*n,(i-1)*n+1:(i)*n) = -Vs{i}';
    
    v((i-2)*n+1:(i-1)*n) = vs{i}'-vs{i-1}';
end

% ATTENTION: DEFINITION OF Gamma IS ONLY FOR n=3!!!!!!!
Gamma = [ 1,  1,  1;
          1,  1, -1;
          1, -1,  1;
          1, -1, -1;
         -1,  1,  1;
         -1,  1, -1;
         -1, -1,  1;
         -1, -1, -1];           

H = kron(eye(N-1),Gamma);
h = ones((N-1)*2^n,1);

Psi = cat(2,invVs{:})';
Vdagger = zeros((N-1)*n, (N-2)*n);

for i = 1:N-2
    Vdagger(1:i*n,(i-1)*n+1:i*n) = Psi(1:i*n,:);
end

% F = [V;
%     -V;
%      H];
%  
% f = [v;
%     -v;
%      h];
%  
% KK = [Vs{1}' zeros(n,(N-2)*n)];
% kk = vs{1}';

gamma = 1.2e9;
kbound = gamma*ones(2^n,1);

% % Nominal LP
% [c,~,flag] = linprog([],H,h,V,v);

% Controller magnitude is below a certain threshold:
[w,~,flag] = linprog(zeros(3,1),[H*Psi; Gamma],[h-H*Vdagger*v;kbound-Gamma*vs{end}'])
% [w,~,flag] = linprog(zeros(3,1),H*Psi,h-H*Vdagger*v);


% [c1,~,flag1] = linprog(ones((N-1)*n,1),H,h,V,v)
% [c2,~,flag2] = linprog(-ones((N-1)*n,1),H,h,V,v)
% [c3,~,flag3] = linprog([],H,h,V,v)

% ccomb = [c1(1:3),c2(1:3),c3(1:3)];
% rank(ccomb)

% KK = blkdiag(Vs{:});
% kk = horzcat(vs{:});
    
if flag==2
    return
end

K = vs{end}+w'
% Gamma*K'

% K = zeros(1,3)


%% Simulations with CT dynamics 

Tsim = 1;

X = zeros(N*n,Tsim/Ts+1);
% X(1:n:end,1) = 47.5+0.5*rand(N,1);
% X(2:n:end,1) = 5+0.666*rand(N,1);
% X(3:n:end,1) = 10+3.9*rand(N,1);
X(1:n:end,1) = 30+30*rand(N,1);
X(2:n:end,1) = 0+9*rand(N,1);
X(3:n:end,1) = 5+10*rand(N,1);

Abig = kron(eye(N),A)-kron(Lp,Ap)+kron(Lc,B*K*1e-0);
eig(Abig);
Bbig = kron(eye(N),B);

expp = @(tau) expm(Abig*(Ts-tau));
integgral = integral(expp,0,Ts,'ArrayValued',true);

time = 0:Ts:Tsim;
for ind = 1:length(time)-1
    X(:,ind+1) = expm(Abig*Ts)*X(:,ind)+integgral*Bbig*48*ones(N,1);
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

subplot(3,1,1)
plot(time,X(1:n:end,:)',line,lw)
grid on
ylabel('$$\mathbf{V_i(V)}$$', font, fn, int,lat,fonts,fs)
legend({'DGU $$1$$','DGU $$2$$','DGU $$3$$','DGU $$4$$','DGU $$5$$','DGU $$6$$','DGU $$7$$','DGU $$8$$','DGU $$9$$'}, font, fn, int,lat,fonts,8,'Orientation','horizontal');
subplot(3,1,2)
plot(time,X(2:n:end,:)',line,lw)
grid on
ylabel('$$\mathbf{I_{ti}(A)}$$', font, fn, int,lat,fonts,fs)
subplot(3,1,3)
plot(time,X(3:n:end,:)',line,lw)
grid on
ylabel('$$\mathbf{v_{i}(Vs)}$$', font, fn, int,lat,fonts,fs)
xlabel('Time ($$\mathbf{s}$$)', font, fn, int,lat,fonts,fs)

% create a new pair of axes inside current figure
axes('position',[.35 .7 .55 .1])
box on % put box around new pair of axes
indexOfInterest = (time < 0.01) & (time >= 0); % range of t near perturbation
plot(time(indexOfInterest),X(1:n:end,indexOfInterest)',line,lw) % plot on new axes
axis tight
% grid on
set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
ylim([30,60])
yticks([30 40 50 60])

% create a new pair of axes inside current figure
axes('position',[.35 .5 .55 .1])
box on % put box around new pair of axes
indexOfInterest = (time < 0.01) & (time >= 0); % range of t near perturbation
plot(time(indexOfInterest),X(2:n:end,indexOfInterest)',line,lw) % plot on new axes
axis tight
% grid on
set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
ylim([-10,25])
yticks([-15 0 15 30])

% create a new pair of axes inside current figure
axes('position',[.35 .3 .55 .1])
box on % put box around new pair of axes
indexOfInterest = (time < 0.01) & (time >= 0); % range of t near perturbation
plot(time(indexOfInterest),X(3:n:end,indexOfInterest)',line,lw) % plot on new axes
axis tight
% grid on
set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
ylim([-5,25])
yticks([0 10 20])





