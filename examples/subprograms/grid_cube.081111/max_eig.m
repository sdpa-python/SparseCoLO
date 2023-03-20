% Frame Optimization with Eigenvalue Constraints
% by using SDP
% 
clear
% 
[dll,matT,nd,nm,coord_x,ir,irr] = member(2);
% 
dmm = zeros(nd,1);
dmm((nd-2):nd,1) = 1.0 * [1; 1; 0];
% dmm((nd-5):(nd-3),1) = 1.0 * [1; 1; 0];
% dmm((nd-8):(nd-6),1) = 1.0 * [1; 1; 0];
% dmm((nd-11):(nd-9),1) = 1.0 * [1; 1; 0];
% dmm((nd-14):(nd-12),1) = 1.0 * [1; 1; 0];
% dmm((nd-17):(nd-15),1) = 1.0 * [1; 1; 0];
% 
cs = (10^(-4)) * (20*30) * ones(nm,1) * 0.01;
cm = cs .* cs;
% 
eeee = 20000;
%%% 200 GPa = 20000 kN/cm^2 = 20000 (10^4 kN)/m^2
rho = 7.86*10^(-4);
%%% 7.86*10^ kg/m^3 = 7.86*10^(-4) (10^4 ton)/m^3
dmm = 2*10^(-3) * dmm;
%%% 10^4 ton = 10^7 kg
% 
%%% make member stiffness matrix : akk(:,i), agg(:,i)
g0  = [1; 0; 0; -1; 0; 0];
d10 = [0; 0; 1; 0; 0; -1];
for i=1:nm
    gg(:,i) = matT(:,:,i)' * g0;
    dd(:,i,1) = matT(:,:,i)' * d10;
    dd(:,i,2) = matT(:,:,i)' *...
        [0;  2*sqrt(3); sqrt(3)*dll(i);...
         0; -2*sqrt(3); sqrt(3)*dll(i)];
end
for i=1:nm
    vec_g(:,i) = sqrt(eeee / dll(i)) * gg(:,i);
    vec_d(:,i,1) = sqrt(eeee / dll(i)) * dd(:,i,1);
    vec_d(:,i,2) = sqrt(eeee / dll(i)^3) * dd(:,i,2);
end
% 
akk = sparse(zeros(nd*nd,nm));
agg = sparse(zeros(nd*nd,nm));
for i=1:nm
    ak = sparse(zeros(nd,nd));
    ak = vec_g(:,i) * vec_g(:,i)';
    akk(:,i) = vec(ak);
end
for i=1:nm
    ag = sparse(zeros(nd,nd));
    for j=1:2
        ag = ag + (vec_d(:,i,j) * vec_d(:,i,j)');
    end
    agg(:,i) = vec(ag) / (4 * pi);
end
%%% <<-- circular cross-section : I = A^(2) / (4 pi)
% 
%%% make member mass matrix : amm(:,i)
%%% Y.W. Kwon & H. bang,
%%% `` The Finite Element Method using MATLAB (2nd ed.)''
%%% pp.282
amm = sparse(zeros(nd*nd,nm));
for i=1:nm
    mat_elm = zeros(6,6);
    % 
    mat_elm(2,2) = 156;
    mat_elm(2,3) = 22 * dll(i);
    mat_elm(2,5) = 54;
    mat_elm(2,6) = -13 * dll(i);
    mat_elm(3,2) = 22 * dll(i);
    mat_elm(3,3) = 4 * (dll(i)^(2));
    mat_elm(3,5) = 13 * dll(i);
    mat_elm(3,6) = -3 * (dll(i)^(2));
    mat_elm(5,2) = 54;
    mat_elm(5,3) = 13 * dll(i);
    mat_elm(5,5) = 156;
    mat_elm(5,6) = -22 * dll(i);
    mat_elm(6,2) = -13 * dll(i);
    mat_elm(6,3) = -3 * (dll(i)^(2));
    mat_elm(6,5) = -22 * dll(i);
    mat_elm(6,6) = 4 * (dll(i)^(2));
    % 
    mat_elm = mat_elm / 420;
    % 
    mat_elm(1,1) = 1/3;
    mat_elm(4,4) = 1/3;
    mat_elm(1,4) = 1/6;
    mat_elm(4,1) = 1/6;
    % 
    mat_elm = (rho * dll(i)) * mat_elm;
    % 
    am = sparse(zeros(nd,nd));
    am = (matT(:,:,i)') * mat_elm * matT(:,:,i);
    amm(:,i) = vec(am);
end
% 
amc = diag(dmm);
% 
fullK = mat((akk * cs) + (agg * cm));
fullM = mat(amm * cs) + amc;
[mode_0,eigV_0] = eig(fullK, fullM);
min_eig = min(diag(eigV_0))

barEig_0 = min_eig;
barV = dll' * cs;

vec_u = 0.05 * mode_0(:,1);

%%%
%%% making SDP data (SeDuMi)...
% ================ for SeDuMi ================
%%% x = [cs, cm, t]
% 
n_dim = (2*nm) + 1;
b = sparse(zeros(n_dim,1));
b(end) = -1;
% 
b = sparse(b);
% 
% ------------ linear inequalities ------------
% 
K.l = [(nm+1)];
% 
%%% nonnegative constraints of areas
% 
At_l_1 = [-eye(nm), zeros(nm,nm+1)];
c_l_1 = zeros(nm,1);
% 
At_l_1 = sparse(At_l_1); c_l_1 = sparse(c_l_1);
% 
%%% volume constraint
% 
At_l_2 = [dll', zeros(1,nm+1)];
c_l_2 = barV;
% 
At_l_2 = sparse(At_l_2);
% 
% -------------- soc inequalities --------------
%%% quadratic inequalities coming from relaxations
% 
At_q_1 = zeros(3*nm,n_dim);
c_q_1 = zeros(3*nm,1);
for i=1:nm
    pp1 = 3 * (i-1);
    pp2 = nm + i;
    At_q_1(pp1+1,pp2) = -1;
    At_q_1(pp1+2,pp2) = -1;
    At_q_1(pp1+3,i) = -1;
    c_q_1(pp1+1,1) = 1/4;
    c_q_1(pp1+2,1) = -1/4;
end
% 
K.q = [3 * ones(1,nm)];
% 
At_q_1 = sparse(At_q_1); c_q_1 = sparse(c_q_1);
% 
%%% quadratic inequality coming from the objective function
% 
At_q_2 = [zeros(1,2*nm), -1;...
        zeros(1,2*nm), -1;...
        -eye(2*nm), zeros(2*nm,1)];
K.q = [K.q, (2*nm+2)];
% 
At_q_2 = sparse(At_q_2);
% 
%%% NOTE: `c_q_2' will be defined in the DCA loop %%%
% 
% -------------- psd inequalities --------------
% 
%%% NOTE: `At_s' & `c_s' will be defined in the DCA loop %%%
% 
K.s = [nd];
% 
%%% parameters...
% ================ for iterations ================
z_ini = [cs; cm];
z_cur = z_ini;
cs_cur = z_cur(1:nm);
cm_cur = z_cur((nm+1):(2*nm));
Vbar_ini = 0.09 * dll' * cs;
Vbar_cur = Vbar_ini;
rho = 0.1;
mu = 0.001;
epsilon = 1.0 * 10^(-4);
his_Eig = [];
his_iter = [];
his_g_star = [];
his_time = [];
z_star = z_cur;
cs_star = cs_cur;
cm_star = cm_cur;
% 
barEig_low = barEig_0;
barEig_upp = barEig_0 * 4.0;
barEig_cur = (barEig_upp + barEig_low)/2;
% 
% % ************  BI-SECTION **** >> comment below >>
while (barEig_upp - barEig_low) > (barEig_0 * (10^(-5)))
% % ************  BI-SECTION **** << comment above <<
    d_norm = 100;
    mu_act = 100;
    i_iter = 0;
    time.DC = cputime;
    while (d_norm > epsilon) & (i_iter < 20)
%     while (i_iter < 10) | (abs(mu_act) > epsilon)
%     while mu_act > epsilon
%     while i_iter < 10
        %%% making SDP data (SeDuMi)...
        % 
        c_q_2 = [(1/4);...
                -(1/4);...
                -(1+(2/rho)) * cs_cur;...
                -cm_cur - (((mu-1)/rho) * ones(nm,1))];
        % 
        At_s = [-akk + (barEig_cur * amm), -agg, zeros(nd*nd,1)];
        c_s = -barEig_cur * vec(amc);
        % 
        At_s = sparse(At_s); c_s = sparse(c_s);
        % 
        % -------------- combine (At, c) --------------
        % 
        At = [At_l_1; At_l_2; At_q_1; At_q_2; At_s];
        c = [c_l_1; c_l_2; c_q_1; c_q_2; c_s];
        % 
        % --------------- execute SeDuMi ---------------
        % 
%         [x,y,info] = sedumi(At, b, c, K);
        pars.fid = 0;
        [x,y,info] = sedumi(At, b, c, K, pars);
        fval = -b' * y;
        z = y(1:(2*nm));
        if (info.pinf + info.dinf) > 0
            break
        end
        % 
        d_norm = norm(z-z_cur,2);
        z_cur = z;
        cs_cur = z_cur(1:nm);
        cm_cur = z_cur((nm+1):(2*nm));
        % 
        cs_m = sqrt(abs(cm_cur));
%         fullK = mat((akk * cs_m) + (agg * cm_cur));
%         fullM = mat(amm * cs_m) + amc;
%         [mode_cur,eigV_cur] = eig(fullK, fullM);
%         min_eig = min(diag(eigV_cur));
        % 
        mu_act = sum((cs_cur.*cs_cur) - cm_cur);
        % 
        gval = (cs_cur' * cs_cur) + ((mu-1) * (ones(1,nm) * cm_cur));
        disp(sprintf('   iter = %g ; mu_act = %g ; d_norm = %g ; gval = %g \n',...
            i_iter, mu_act, d_norm, gval));
        % 
        i_iter = i_iter + 1;
    end
    time.DC = cputime - time.DC;
    his_Eig = [his_Eig, barEig_cur];
    his_iter = [his_iter, i_iter];
    his_g_star = [his_g_star, gval];
    his_time = [his_time, time.DC];
    % 
    g_star = gval;
    % 
% % ************  BI-SECTION **** >> comment below >>
    if g_star > 0
        barEig_low = barEig_cur;
        barEig_cur = (barEig_upp + barEig_cur)/2;
        % 
        z_star = z_cur;
        cs_star = cs_cur;
        cm_star = cm_cur;
    else
        barEig_upp = barEig_cur;
        barEig_cur = (barEig_low + barEig_cur)/2;
    end
    disp(sprintf(' ###### \n\n'));
%     [barEig_upp, barEig_low]
    % 
    z_cur = z_star;
    cs_cur = cs_star;
    cm_cur = cm_star;
end
% % ************  BI-SECTION **** << comment above <<


barEig = barEig_cur

fullK = mat((akk * cs_cur) + (agg * cm_cur));
fullM = mat(amm * cs_cur) + amc;
[mode_cur,eigV_cur] = eig(fullK, fullM);
eigV_cur = diag(eigV_cur);
eigV_cur(1:5)'

% mu_act = scale0 * sum((cs_cur.*cs_cur) - cm_cur)

fullK = mat((akk * cs_cur) + (agg * (cs_cur .* cs_cur)));
fullM = mat(amm * cs_cur) + amc;
[mode_cs,eigV_cs] = eig(fullK, fullM);
eigV_cs = diag(eigV_cs);
eig_tmp = eigV_cs - (barEig * ones(nd,1));
[min_eig_cs, ind_eig] = min(abs(eig_tmp));
min_eig_cs = eigV_cs(ind_eig:(ind_eig+4))'

cs0 = cs_cur;
for i=1:nm
    if cs0(i) < 10^(-4)
        cs0(i) = 10^(-12);
    end
end
fullK0 = mat((akk * cs0) + (agg * (cs0 .* cs0)));
fullM0 = mat(amm * cs0) + amc;
[mode_cs0,eigV_cs0] = eig(fullK0, fullM0);
eigV_cs0 = diag(eigV_cs0);
eig_tmp = eigV_cs0 - (10^(-3) * ones(nd,1));
ind_eig0 = 1;
for i=1:nd
    if eig_tmp(i) < 0
        ind_eig0 = ind_eig0 + 1;
    end
end
min_eig_cs0 = eigV_cs0((ind_eig0):(ind_eig0+6))'

vec_u = 0.0001 * mode_cur(:,3);
% vec_u = 0.002 * mode_cs(:,ind_eig+1);
% dummy = draw_mode(vec_u,coord_x,irr,ir,dll,nm,matT,cs_cur);

vec_u = 0.01 * mode_cs0(:,ind_eig0);
dummy = draw_mode(vec_u,coord_x,irr,ir,dll,nm,matT,cs0);

% [cs_cur, cs_m, cm_cur]

%%%% save optimal cross-sections
qqq = [99; cs_cur];
save cs0.dat -ascii qqq


