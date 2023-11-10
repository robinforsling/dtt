function [COIN,ANEES,GEVO,PCO] = gevo_evaluation_fusion_example
% --- gevo_evaluation_fusion_example() ------------------------------------
% Section 5.3.7 Parametrized Fusion Example - GEVO: Figure 5.9
%
% 2023-10-30 Robin Forsling

set_latex_interpreter;

case_common_info = 4;


% --- PARAMETERS ---
nx = 6;
Ai = diag([64 32 16 8 4 2]);
Bi = diag([5 8 13 21 34 55]);
Ci = 2*get_common_information(case_common_info);

m_vec = [1 2 3];
nm = length(m_vec);

rho_vec = 0.025:0.05:0.975; %0.01:0.02:0.99;
nrho = length(rho_vec);


% --- DATA RECORDS ---
bsc_gevo = zeros(nm,nrho); 
dkf_gevo = bsc_gevo; ci_gevo = bsc_gevo; le_gevo = bsc_gevo; nkf_gevo = bsc_gevo;
dkf_pco = bsc_gevo; ci_pco = bsc_gevo; le_pco = bsc_gevo; nkf_pco = bsc_gevo;
s = zeros(nm,nrho); dkf_c.anees = s; dkf_c.coin = s; 
ci_c = dkf_c; le_c = dkf_c; nkf_c = dkf_c;


% --- RUN SIMULATIONS ---
for i = 1:nm

    m = m_vec(i);

    for j = 1:nrho
    
        rho = rho_vec(j);   

        Gi = rho*Ci;
        R1 = inv((1-rho)*Ai + Gi); 
        R2 = inv((1-rho)*Bi + Gi);
        R2d = inv((1-rho)*Bi);
        R12 = R1*Gi*R2; 
        R12d = zeros(nx);

           
        % GEVO
        Psi_bsc = gevobsc_dr(R1,R2,R12,m);
        Psi_dkf = gevokf_dr(R1,R2d,m);
        Psi_ci = gevoci_dr(R1,R2,m);
        Psi_le = gevole_dr(R1,R2,m);
        Psi_nkf = gevokf_dr(R1,R2,m);
    
        Pbsc = bsc_true(R1,R2,R12,Psi_bsc);
        [Pdkf,Pdkf0] = kf_true(R1,R2d,R12d,Psi_dkf);
        [Pci,Pci0] = ci_true(R1,R2,R12,Psi_ci);
        [Ple,Ple0] = le_true(R1,R2,R12,Psi_le);
        [Pnkf,Pnkf0] = kf_true(R1,R2,R12,Psi_nkf);

        bsc_gevo(i,j) = sqrt(trace(Pbsc));
        dkf_gevo(i,j) = sqrt(trace(Pdkf));
        ci_gevo(i,j) = sqrt(trace(Pci));
        le_gevo(i,j) = sqrt(trace(Ple));
        nkf_gevo(i,j) = sqrt(trace(Pnkf));

        c_dkf = check_conservativeness(Pdkf,Pdkf0);
        c_ci = check_conservativeness(Pci,Pci0);
        c_le = check_conservativeness(Ple,Ple0);
        c_nkf = check_conservativeness(Pnkf,Pnkf0);

        dkf_c.anees(i,j) = c_dkf.anees; dkf_c.coin(i,j) = c_dkf.coin; 
        ci_c.anees(i,j) = c_ci.anees; ci_c.coin(i,j) = c_ci.coin;
        le_c.anees(i,j) = c_le.anees; le_c.coin(i,j) = c_le.coin;
        nkf_c.anees(i,j) = c_nkf.anees; nkf_c.coin(i,j) = c_nkf.coin;


        % PCO
        Psi_pco_d = pco_dr(R2d,m);
        Psi_pco = pco_dr(R2,m);

        [Pdkf,Pdkf0] = kf_true(R1,R2d,R12d,Psi_pco_d);
        [Pci,Pci0] = ci_true(R1,R2,R12,Psi_pco);
        [Ple,Ple0] = le_true(R1,R2,R12,Psi_pco);
        [Pnkf,Pnkf0] = kf_true(R1,R2,R12,Psi_pco);

        dkf_pco(i,j) = sqrt(trace(Pdkf));
        ci_pco(i,j) = sqrt(trace(Pci));
        le_pco(i,j) = sqrt(trace(Ple));
        nkf_pco(i,j) = sqrt(trace(Pnkf));
    end
end


% --- STORE DATA ---
COIN.dkf = dkf_c.coin; 
COIN.ci = ci_c.coin;
COIN.le = le_c.coin;
COIN.nkf = nkf_c.coin;

ANEES.dkf = dkf_c.anees; 
ANEES.ci = ci_c.anees;
ANEES.le = le_c.anees;
ANEES.nkf = nkf_c.anees;

GEVO.bsc = bsc_gevo;
GEVO.dkf = dkf_gevo;
GEVO.ci = ci_gevo;
GEVO.le = le_gevo;
GEVO.nkf = nkf_gevo;

PCO.dkf = dkf_pco;
PCO.ci = ci_pco;
PCO.le = le_pco;
PCO.nkf = nkf_pco;


% --- RESULTS ---
clr = get_thesis_colors;
clrdkf = clr.darkyellow;

lw_vec = 1.0:0.5:4;
ls = {'-','--',':','-.','-','--','-.',':'};

figure(1);clf


% --- COIN --
subplot(2,2,1);hold on

for i = 1:nm; h = plot(rho_vec,dkf_c.coin(i,:),'-','DisplayName',sprintf('dKF $m=%d$',m_vec(i))); set_h(h,clrdkf,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,ci_c.coin(i,:),'-','DisplayName',sprintf('CI $m=%d$',m_vec(i))); set_h(h,clr.ci,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,le_c.coin(i,:),'-','DisplayName',sprintf('LE $m=%d$',m_vec(i))); set_h(h,clr.le,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,nkf_c.coin(i,:),'-','DisplayName',sprintf('nKF $m=%d$',m_vec(i))); set_h(h,clr.naive,lw_vec(i),ls{i}); end

xlim([0 1]); ylim([0.5 2.0])
legend('show'); 
box on; grid on
xlabel('$\rho$','interpreter','latex'); 
ylabel('COIN','interpreter','latex'); 


% --- ANEES ---    
subplot(2,2,2);hold on

for i = 1:nm; h = plot(rho_vec,dkf_c.anees(i,:),'-','DisplayName',sprintf('dKF $m=%d$',m_vec(i))); set_h(h,clrdkf,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,ci_c.anees(i,:),'-','DisplayName',sprintf('CI $m=%d$',m_vec(i))); set_h(h,clr.ci,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,le_c.anees(i,:),'-','DisplayName',sprintf('LE $m=%d$',m_vec(i))); set_h(h,clr.le,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,nkf_c.anees(i,:),'-','DisplayName',sprintf('nKF $m=%d$',m_vec(i))); set_h(h,clr.naive,lw_vec(i),ls{i}); end

xlim([0 1]); ylim([0.5 2.0])
box on; grid on
xlabel('$\rho$','interpreter','latex'); 
ylabel('ANEES','interpreter','latex'); 


% --- RMTR GEVO-TO-BSC ---
subplot(2,2,3);hold on
 
for i = 1:nm; h = plot(rho_vec,dkf_gevo(i,:)./bsc_gevo(i,:),'-','DisplayName',sprintf('dKF $m=%d$',m_vec(i))); set_h(h,clrdkf,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,ci_gevo(i,:)./bsc_gevo(i,:),'-','DisplayName',sprintf('CI $m=%d$',m_vec(i))); set_h(h,clr.ci,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,le_gevo(i,:)./bsc_gevo(i,:),'-','DisplayName',sprintf('LE $m=%d$',m_vec(i))); set_h(h,clr.le,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,nkf_gevo(i,:)./bsc_gevo(i,:),'-','DisplayName',sprintf('nKF $m=%d$',m_vec(i))); set_h(h,clr.naive,lw_vec(i),ls{i}); end

xlim([0 1]); ylim([0.75 1.5])
box on; grid on
xlabel('$\rho$','interpreter','latex'); 
ylabel('RMTR-GTB','interpreter','latex'); 


% --- RMTR GEVO-TO-PCO ---
subplot(2,2,4);hold on
 
for i = 1:nm; h = plot(rho_vec,dkf_gevo(i,:)./dkf_pco(i,:),'-','DisplayName',sprintf('dKF $m=%d$',m_vec(i))); set_h(h,clrdkf,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,ci_gevo(i,:)./ci_pco(i,:),'-','DisplayName',sprintf('CI $m=%d$',m_vec(i))); set_h(h,clr.ci,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,le_gevo(i,:)./le_pco(i,:),'-','DisplayName',sprintf('LE $m=%d$',m_vec(i))); set_h(h,clr.le,lw_vec(i),ls{i}); end
for i = 1:nm; h = plot(rho_vec,nkf_gevo(i,:)./nkf_pco(i,:),'-','DisplayName',sprintf('nKF $m=%d$',m_vec(i))); set_h(h,clr.naive,lw_vec(i),ls{i}); end

xlim([0 1]); ylim([0.775 1.025])
box on; grid on
xlabel('$\rho$','interpreter','latex'); 
ylabel('RMTR-GTP','interpreter','latex');

set_fontsize_all(14)


end




% --- GEVO & PCO ----------------------------------------------------------
function Psi = gevobsc_dr(R1,R2,R12,m)
    Q = make_symmetric((R1-R12)'*(R1-R12)); 
    S = make_symmetric(R1+R2-R12-R12');
    [X,G] = eig(Q,S,'qz');
    idx = get_max_idx_vec(G,m);
    Psi = gram_schmidt_process(X(:,idx)');
end

function Psi = gevokf_dr(R1,R2,m)
    Q = R1*R1;
    S = R1+R2;
    [X,G] = eig(Q,S);
    idx = get_max_idx_vec(G,m);
    Psi = gram_schmidt_procedure(X(:,idx)');
end

function Psi = gevoci_dr(R1,R2,m)
    max_iter = 10; epsilon = 1e-4;
    wmin = 0.001; wmax = 1;
    I1 = inv(R1); 
    w = 0.5; J = Inf;
    
    for k = 1:max_iter
    
        J_prev = J;
    
        A = R1/w; B = R2/(1-w);
        Q = A^2; S = A+B;
        [X,G] = eig(Q,S);
        idx = get_max_idx_vec(G,m);
        Psi = X(:,idx)';
        RPsi = Psi*R2*Psi';
        
        f = @(w) trace( inv(w*I1 + (1-w)*Psi'/RPsi*Psi) );
        w = fminbnd(f,wmin,wmax);
        J = gevoci_loss_function(I1,RPsi,Psi,w);
    
        if abs(J_prev-J)/J < epsilon; break; end
    end
end

function Psi = gevole_dr(R1,R2,m)
    nx = size(R1,1); 
    [U1,D1] = eig(make_symmetric(inv(R1))); 
    T1 = inv(sqrtm(D1))*U1';
    [U2,D2] = eig(make_symmetric(T1/R2*T1'));
    T = U2'*T1; Ti = inv(T);
    Gammai = Ti*diag(min([ones(1,nx) ; diag(D2)']))*Ti';
    R12 = R1*Gammai*R2; R12 = 0.99*R12;
    Q = make_symmetric((R1-R12)'*(R1-R12)); 
    S = make_symmetric(R1+R2-R12-R12');
    [X,G] = eig(Q,S,'qz');
    idx = get_max_idx_vec(G,m);
    Psi = gram_schmidt_process(X(:,idx)');
end

function J = gevoci_loss_function(I1,RPsi,Psi,w)
    J = trace( inv(w*I1+(1-w)*Psi'/RPsi*Psi) );
end

function Psi = pco_dr(R2,m)
    [U,D] = eig(R2);
    idx = get_min_idx_vec(D,m);
    Psi = U(:,idx)';
end

function U = gram_schmidt_procedure(X)
    [n,m] = size(X);
    is_tall = n > m;    
    if ~is_tall; X = X'; [n,m] = size(X); end 
    
    U = zeros(n,m);
    U(:,1) = X(:,1)/norm(X(:,1));
    
    for i = 2:m
        x = X(:,i); u = x;
        for j = 1:i-1
            uprev = U(:,j);
            u = u - po(uprev,x);
        end
        U(:,i) = u/norm(u);
    end    
    if ~is_tall; U = U'; end  
end

function vu = po(u,v)
    vu = (u'*v/(u'*u))*u;
end




% --- EVALUATION ----------------------------------------------------------
function c = check_conservativeness(P,P0)
    nx = size(P,1);
    c.anees = trace(P0/P)/nx;
    L = chol(P,'lower'); Li = inv(L);
    c.coin = max(eig(Li*P0*Li'));
end

function P = bsc_true(R1,R2,R12,Psi)
    D = R1-R12; S = R1+R2-R12-R12';
    P = R1 - D*Psi'/(Psi*S*Psi')*Psi*D';
end

function [P,P0] = kf_true(R1,R2,R12,Psi)
    RPsi = Psi*R2*Psi';
    R12Psi = R12*Psi';
    R = [R1 R12Psi ; R12Psi' RPsi];
    P = inv(inv(R1) + Psi'/RPsi*Psi);
    K1 = P/R1; K2 = P*Psi'/RPsi;
    K = [K1 K2];
    P0 = K*R*K';
end

function [P,P0] = ci_true(R1,R2,R12,Psi)
    RPsi = Psi*R2*Psi';
    R12Psi = R12*Psi';   
    f = @(w) trace(inv(w*inv(R1) + (1-w)*Psi'/RPsi*Psi));
    w = fminbnd(f,0,1);
    P = inv(w*inv(R1) + (1-w)*Psi'/RPsi*Psi);
    R0 = [R1 R12Psi ; R12Psi' RPsi]; 
    K1 = w*P/R1; K2 = (1-w)*P*Psi'/RPsi;
    K = [K1 K2];
    P0 = K*R0*K';
end

function [P,P0] = le_true(R1,R2,R12,Psi)
    nx = size(R1,1);  
    RPsi = Psi*R2*Psi'; R12Psi = R12*Psi';
    R = [R1 R12Psi ; R12Psi' RPsi];
    I1 = inv(R1); IPsi = inv(RPsi); I2 = Psi'*IPsi*Psi;
    [U1,D1] = eig(make_symmetric(I1)); T1 = inv(sqrtm(D1))*U1';
    [U2,D2] = eig(make_symmetric(T1*I2*T1')); T = U2'*T1; Ti = inv(T);
    P = inv(Ti*diag(max([ones(1,nx) ; diag(D2)']))*Ti');
    Gammai = Ti*diag(min([ones(1,nx) ; diag(D2)']))*Ti';
    R12 = R1*Gammai*R2;
    S = R1+R2-R12-R12';
    K2 = (R1-R12)*Psi'/(Psi*S*Psi'); 
    K1 = eye(nx) - K2*Psi;
    K = [K1 K2];
    P0 = K*R*K';
end



% --- OTHER STUFF ---
function Ci = get_common_information(idx)
    switch idx
        case 1; Ci = [10 2 2 1 -1 1 ; 2 10 4 -4 -2 -2 ; 2 4 15 0 -2 -2 ; 1 -4 0 15 -1 0 ; -2 -2 -2 -1 5 0 ; 1 -2 -2 0 0 5];
        case 2; Ci = [10 2 2 0 0 0 ; 2 10 4 -4 -2 -2 ; 2 4 15 0 -2 -2 ; 0 -4 0 15 0 0 ; -2 -2 -2 0 5 0 ; 0 -2 -2 0 0 5];
        case 3; Ci = [8 2 2 1 -1 1 ; 2 10 4 -4 -2 -2 ; 2 4 15 0 -2 -2 ; 1 -4 0 25 -1 0 ; -2 -2 -2 -1 5 0 ; 1 -2 -2 0 0 10];
        case 4; Ci = [8 2 2 0 0 0 ; 2 10 4 -4 -2 -2 ; 2 4 15 0 -2 -2 ; 0 -4 0 25 0 0 ; -2 -2 -2 0 5 0 ; 0 -2 -2 0 0 10]; 
    end
    Ci = make_symmetric(Ci);
end

function set_h(h,clr,lw,ls)
    h.Color = clr; 
    h.LineWidth = lw; 
    h.LineStyle = ls;
end



