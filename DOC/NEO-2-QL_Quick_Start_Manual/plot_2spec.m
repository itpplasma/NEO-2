%% define constants
e=4.8032e-10; % elementary charge

%% load NEO-2 output
% 2 species (ExB only)
data_NEO2_nspec2_lag7_newint_wHelCore=h52struct('final_neo2_multispecies_out.h5');
fname_NEO2_nspec2_lag7_newint_wHelCore=fieldnames(data_NEO2_nspec2_lag7_newint_wHelCore);

% allocate storage array
num_data_NEO2=1; % total number of HDFF5 data files to be processed
data_NEO2=cell(num_data_NEO2,2);
% 2 species (ExB only)
data_NEO2{1,1}=data_NEO2_nspec2_lag7_newint_wHelCore;
data_NEO2{1,2}=fname_NEO2_nspec2_lag7_newint_wHelCore;

% extract data
mspec_NEO2=cell(num_data_NEO2,1);
zspec_NEO2=cell(num_data_NEO2,1);
boozer_s_NEO2=cell(num_data_NEO2,1);
boozer_psi_pr_NEO2=cell(num_data_NEO2,1);
aiota_NEO2=cell(num_data_NEO2,1);
bcovar_tht_NEO2=cell(num_data_NEO2,1);
bcovar_phi_NEO2=cell(num_data_NEO2,1);
avb2_NEO2=cell(num_data_NEO2,1);
Er_NEO2=cell(num_data_NEO2,1);
avEparB_ov_avb2_NEO2=cell(num_data_NEO2,1);
avnabstor_NEO2=cell(num_data_NEO2,1);
Mt_e_NEO2=cell(num_data_NEO2,1);
Mt_d_NEO2=cell(num_data_NEO2,1);
TphiNA_tot_NEO2=cell(num_data_NEO2,1);
TphiNA_spec_NEO2=cell(num_data_NEO2,1);
Gamma_AX_spec_NEO2=cell(num_data_NEO2,1);
Gamma_AX_Ware_spec_NEO2=cell(num_data_NEO2,1);
Gamma_NA_spec_NEO2=cell(num_data_NEO2,1);
Gamma_NA_Ware_spec_NEO2=cell(num_data_NEO2,1);
n_spec_NEO2=cell(num_data_NEO2,1);
T_spec_NEO2=cell(num_data_NEO2,1);
Qflux_AX_spec_NEO2=cell(num_data_NEO2,1);
Qflux_AX_Ware_spec_NEO2=cell(num_data_NEO2,1);
Qflux_NA_spec_NEO2=cell(num_data_NEO2,1);
Qflux_NA_Ware_spec_NEO2=cell(num_data_NEO2,1);
for file_ind=1:num_data_NEO2
    % get structure
    data_struct=data_NEO2{file_ind,1};
    data_fname=data_NEO2{file_ind,2};
    % get NEO-2 data
    data_ctr=0;
    mspec_NEO2_tmp=cell(numel(data_fname),1);
    zspec_NEO2_tmp=cell(numel(data_fname),1);
    boozer_s_NEO2_tmp=zeros(numel(data_fname),1);
    psi_pr_hat_NEO2_tmp=zeros(numel(data_fname),1);
    Bref_NEO2_tmp=zeros(numel(data_fname),1);
    aiota_NEO2_tmp=zeros(numel(data_fname),1);
    bcovar_tht_NEO2_tmp=zeros(numel(data_fname),1);
    bcovar_phi_NEO2_tmp=zeros(numel(data_fname),1);
    avbhat2_NEO2_tmp=zeros(numel(data_fname),1);
    R0_NEO2_tmp=zeros(numel(data_fname),1);
    Er_NEO2_tmp=zeros(numel(data_fname),1);
    avEparB_ov_avb2_NEO2_tmp=zeros(numel(data_fname),1);
    avnabstor_NEO2_tmp=zeros(numel(data_fname),1);
    MteOvR_NEO2_tmp=zeros(numel(data_fname),1);
    MtdOvR_NEO2_tmp=zeros(numel(data_fname),1);
    TphiNA_tot_NEO2_tmp=zeros(numel(data_fname),1);
    TphiNA_spec_NEO2_tmp=cell(numel(data_fname),1);
    Gamma_AX_spec_NEO2_tmp=cell(numel(data_fname),1);
    Gamma_AX_Ware_spec_NEO2_tmp=cell(numel(data_fname),1);
    Gamma_NA_spec_NEO2_tmp=cell(numel(data_fname),1);
    Gamma_NA_Ware_spec_NEO2_tmp=cell(numel(data_fname),1);
    n_spec_NEO2_tmp=cell(numel(data_fname),1);
    T_spec_NEO2_tmp=cell(numel(data_fname),1);
    Qflux_AX_spec_NEO2_tmp=cell(numel(data_fname),1);
    Qflux_AX_Ware_spec_NEO2_tmp=cell(numel(data_fname),1);
    Qflux_NA_spec_NEO2_tmp=cell(numel(data_fname),1);
    Qflux_NA_Ware_spec_NEO2_tmp=cell(numel(data_fname),1);
    for es_ind=1:numel(data_fname)
        if (~strcmp(data_fname{es_ind}(1:3),'es_'))
            continue
        end
        data_ctr=data_ctr+1;
        mspec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).m_spec;
        zspec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).z_spec;
        es_str=data_fname{es_ind}(4:end);
        boozer_s_NEO2_tmp(es_ind)=str2num(strrep(es_str,'p','.'));
        psi_pr_hat_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).psi_pr_hat;
        Bref_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).Bref;
        aiota_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).aiota;
        bcovar_tht_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).bcovar_tht;
        bcovar_phi_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).bcovar_phi;
        avbhat2_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).avbhat2;
        R0_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).R0;
        Er_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).Er;
        avEparB_ov_avb2_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).avEparB_ov_avb2;
        avnabstor_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).avnabstor;
        MteOvR_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).MtOvR(1);
        MtdOvR_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).MtOvR(2);
        TphiNA_tot_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).TphiNA_tot;
        TphiNA_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).TphiNA_spec;
        Gamma_AX_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Gamma_AX_spec;
        Gamma_AX_Ware_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Gamma_AX_Ware_spec;
        Gamma_NA_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Gamma_NA_spec;
        Gamma_NA_Ware_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Gamma_NA_Ware_spec;
        n_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).n_spec;
        T_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).T_spec;
        Qflux_AX_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Qflux_AX_spec;
        Qflux_AX_Ware_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Qflux_AX_Ware_spec;
        Qflux_NA_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Qflux_NA_spec;
        Qflux_NA_Ware_spec_NEO2_tmp{es_ind}=data_struct.(data_fname{es_ind}).Qflux_NA_Ware_spec;
    end
    % store data for plotting
    mspec_NEO2{file_ind}=mspec_NEO2_tmp(1:data_ctr);
    zspec_NEO2{file_ind}=zspec_NEO2_tmp(1:data_ctr);
    boozer_s_NEO2{file_ind}=boozer_s_NEO2_tmp(1:data_ctr);
    boozer_psi_pr_NEO2{file_ind}=psi_pr_hat_NEO2_tmp(1:data_ctr).*Bref_NEO2_tmp(1:data_ctr);
    aiota_NEO2{file_ind}=aiota_NEO2_tmp(1:data_ctr);
    bcovar_tht_NEO2{file_ind}=bcovar_tht_NEO2_tmp(1:data_ctr);
    bcovar_phi_NEO2{file_ind}=bcovar_phi_NEO2_tmp(1:data_ctr);
    avb2_NEO2{file_ind}=avbhat2_NEO2_tmp(1:data_ctr).*(Bref_NEO2_tmp(1:data_ctr).^2);
    Er_NEO2{file_ind}=Er_NEO2_tmp(1:data_ctr);
    avEparB_ov_avb2_NEO2{file_ind}=avEparB_ov_avb2_NEO2_tmp(1:data_ctr);
    avnabstor_NEO2{file_ind}=avnabstor_NEO2_tmp(1:data_ctr);
    Mt_e_NEO2{file_ind}=MteOvR_NEO2_tmp(1:data_ctr).*R0_NEO2_tmp(1:data_ctr);
    Mt_d_NEO2{file_ind}=MtdOvR_NEO2_tmp(1:data_ctr).*R0_NEO2_tmp(1:data_ctr);
    TphiNA_tot_NEO2{file_ind}=TphiNA_tot_NEO2_tmp(1:data_ctr);
    TphiNA_spec_NEO2{file_ind}=TphiNA_spec_NEO2_tmp(1:data_ctr);
    Gamma_AX_spec_NEO2{file_ind}=Gamma_AX_spec_NEO2_tmp(1:data_ctr);
    Gamma_AX_Ware_spec_NEO2{file_ind}=Gamma_AX_Ware_spec_NEO2_tmp(1:data_ctr);
    Gamma_NA_spec_NEO2{file_ind}=Gamma_NA_spec_NEO2_tmp(1:data_ctr);
    Gamma_NA_Ware_spec_NEO2{file_ind}=Gamma_NA_Ware_spec_NEO2_tmp(1:data_ctr);
    n_spec_NEO2{file_ind}=n_spec_NEO2_tmp(1:data_ctr);
    T_spec_NEO2{file_ind}=T_spec_NEO2_tmp(1:data_ctr);
    Qflux_AX_spec_NEO2{file_ind}=Qflux_AX_spec_NEO2_tmp(1:data_ctr);
    Qflux_AX_Ware_spec_NEO2{file_ind}=Qflux_AX_Ware_spec_NEO2_tmp(1:data_ctr);
    Qflux_NA_spec_NEO2{file_ind}=Qflux_NA_spec_NEO2_tmp(1:data_ctr);
    Qflux_NA_Ware_spec_NEO2{file_ind}=Qflux_NA_Ware_spec_NEO2_tmp(1:data_ctr);
end

%% Compute integral NTV torque and surface area

TphiNA_int_tot_NEO2=cell(num_data_NEO2,1);
surf_area_NEO2=cell(num_data_NEO2,1);
for file_ind=1:num_data_NEO2
    % input
    boozer_s_NEO2_tmp=boozer_s_NEO2{file_ind};
    boozer_psi_pr_NEO2_tmp=boozer_psi_pr_NEO2{file_ind};
    aiota_NEO2_tmp=aiota_NEO2{file_ind};
    bcovar_tht_NEO2_tmp=bcovar_tht_NEO2{file_ind};
    bcovar_phi_NEO2_tmp=bcovar_phi_NEO2{file_ind};
    avb2_NEO2_tmp=avb2_NEO2{file_ind};
    avnabstor_NEO2_tmp=avnabstor_NEO2{file_ind};
    TphiNA_tot_NEO2_tmp=TphiNA_tot_NEO2{file_ind};
    % local array
    TphiNA_int_tot_NEO2_tmp=zeros(size(TphiNA_tot_NEO2_tmp));
    % compute integral torque
    TphiNA_int_tot_NEO2_tmp(1)=0;
    for es_ind=2:numel(TphiNA_tot_NEO2_tmp)
        TphiNA_int_tot_NEO2_tmp(es_ind)=...
            trapz(boozer_s_NEO2_tmp(1:es_ind),TphiNA_tot_NEO2_tmp(1:es_ind).*...
            (aiota_NEO2_tmp(1:es_ind).*bcovar_tht_NEO2_tmp(1:es_ind)+bcovar_phi_NEO2_tmp(1:es_ind))./...
            avb2_NEO2_tmp(1:es_ind));
    end
    TphiNA_int_tot_NEO2_tmp=(4*pi^2)*boozer_psi_pr_NEO2_tmp.*TphiNA_int_tot_NEO2_tmp;
    % store data for plotting
    TphiNA_int_tot_NEO2{file_ind}=TphiNA_int_tot_NEO2_tmp;
    surf_area_NEO2{file_ind}=(4*(pi^2)*avnabstor_NEO2_tmp./avb2_NEO2_tmp).*...
        abs(boozer_psi_pr_NEO2_tmp.*(aiota_NEO2_tmp.*bcovar_tht_NEO2_tmp+bcovar_phi_NEO2_tmp));
end

%% Compute radial and inductive electric field for Hakan Smith

% 2 species results (3D VMEC equilibrium 31021 + Profiles 32169 wo hel. core)
ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

E_rhotor=2*Er_NEO2{ind_data_NEO2}.*sqrt(boozer_s_NEO2{ind_data_NEO2})./avnabstor_NEO2{ind_data_NEO2}; % [statV]
E_rhotor=E_rhotor*3e2*1e-3; % [kV]

E_rhotor_for_Hakan=[sqrt(boozer_s_NEO2{ind_data_NEO2}),E_rhotor];

avEparB_ov_SqrtAvB2=avEparB_ov_avb2_NEO2{ind_data_NEO2}.*sqrt(avb2_NEO2{ind_data_NEO2});
avEparB_ov_SqrtAvB2=avEparB_ov_SqrtAvB2*3e4*1e3; % [mV/m]

EparB_for_Hakan=[sqrt(boozer_s_NEO2{ind_data_NEO2}),avEparB_ov_SqrtAvB2]

figure(2)
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),avEparB_ov_SqrtAvB2)
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\langle E_{\parallel} B \rangle \ / \ \langle B^2 \rangle^{1/2}$ [mV/m]',...
    'interpreter','latex')
set(gca,'XTick',[0:0.2:1])
export_plot('avEparB_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_Spec2')

%return
%% 2 species results (3D VMEC equilibrium 32138 + Profiles 32169 w/ hel. core) 

ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

% surface area for computation of total fluxes
surf_area=surf_area_NEO2{ind_data_NEO2};

% Get NEO-2 data for plots
num_radial_pts=numel(Gamma_AX_spec_NEO2{ind_data_NEO2});
Gamma_AX_e=zeros(num_radial_pts,1);
Gamma_AX_d=zeros(num_radial_pts,1);
Gamma_AX_Ware_e=zeros(num_radial_pts,1);
Gamma_AX_Ware_d=zeros(num_radial_pts,1);
Qflux_AX_e=zeros(num_radial_pts,1);
Qflux_AX_d=zeros(num_radial_pts,1);
Qflux_AX_Ware_e=zeros(num_radial_pts,1);
Qflux_AX_Ware_d=zeros(num_radial_pts,1);
ne=zeros(num_radial_pts,1);
nd=zeros(num_radial_pts,1);
Te=zeros(num_radial_pts,1);
Td=zeros(num_radial_pts,1);
j_AX_e=zeros(num_radial_pts,1);
j_AX_d=zeros(num_radial_pts,1);
j_AX_Ware_e=zeros(num_radial_pts,1);
j_AX_Ware_d=zeros(num_radial_pts,1);
for k=1:num_radial_pts
    Gamma_AX_e(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_AX_d(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_AX_Ware_e(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_AX_Ware_d(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_AX_e(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_AX_d(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_AX_Ware_e(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_AX_Ware_d(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    ne(k)=n_spec_NEO2{ind_data_NEO2}{k}(1);
    nd(k)=n_spec_NEO2{ind_data_NEO2}{k}(2);
    Te(k)=T_spec_NEO2{ind_data_NEO2}{k}(1);
    Td(k)=T_spec_NEO2{ind_data_NEO2}{k}(2);
    j_AX_e(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_AX_d(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_AX_Ware_e(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_AX_Ware_d(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
end

Qflux_AX_e_tot=Qflux_AX_e.*surf_area;
Qflux_AX_d_tot=Qflux_AX_d.*surf_area;

j_AX_e_tot=j_AX_e.*surf_area;
j_AX_d_tot=j_AX_d.*surf_area;

Gamma_AX_woWare_e=Gamma_AX_e-Gamma_AX_Ware_e;
Gamma_AX_woWare_d=Gamma_AX_d-Gamma_AX_Ware_d;

Qflux_AX_woWare_e=Qflux_AX_e-Qflux_AX_Ware_e;
Qflux_AX_woWare_d=Qflux_AX_d-Qflux_AX_Ware_d;
Qflux_AX_woWare_e_tot=Qflux_AX_woWare_e.*surf_area;
Qflux_AX_woWare_d_tot=Qflux_AX_woWare_d.*surf_area;

j_AX_woWare_e=j_AX_e-j_AX_Ware_e;
j_AX_woWare_d=j_AX_d-j_AX_Ware_d;
j_AX_woWare_e_tot=j_AX_woWare_e.*surf_area;
j_AX_woWare_d_tot=j_AX_woWare_d.*surf_area;

figure(1000)
lb_1000=0;
rb_1000=0.9;
tb_1000=0.03;
bb_1000=0.0;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_d./nd),'r:>')  % conversion SI units
%plot([lb_1000,rb_1000],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1000,rb_1000])
ylim([bb_1000,tb_1000])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
legend('NEO-2: e ','NEO-2: d','Location','North')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Gamma_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_AllCodes')

figure(1001)
lb_1001=0;
rb_1001=0.9;
tb_1001=0.03;
bb_1001=-0.06;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_e./ne),'b:o') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_d./nd),'r:+')  % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_e./ne),'b:sq') % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_d./nd),'r:>')  % conversion SI units
plot([lb_1001,rb_1001],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1001,rb_1001])
ylim([bb_1001,tb_1001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
lh=legend('e: w/ $E_\parallel$','d: w/ $E_\parallel$',...
    'e: w/o $E_\parallel$','d: w/o $E_\parallel$','Location','SouthWest');
set(lh,'interpreter','latex')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Gamma_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_NEO2')

figure(1010)
lb_1010=0;
ub_1010=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_woWare_d./nd),'r:>')  % conversion SI units
%plot([lb_1010,ub_1010],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1010,ub_1010])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_AllCodes')

figure(10100)
lb_10100=0;
ub_10100=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_AX_woWare_e_tot,'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_AX_woWare_d_tot,'r:>')  % conversion SI units
%plot([lb_10100,ub_10100],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10100,ub_10100])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_tot_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_AllCodes')

figure(1011)
lb_1011=0;
ub_1011=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_d./nd),'r:>')  % conversion SI units
%plot([lb_1011,ub_1011],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1011,ub_1011])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_AllCodes')

figure(10110)
lb_10110=0;
ub_10110=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_AX_e_tot,'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_AX_d_tot,'r:>')  % conversion SI units
%plot([lb_10110,ub_10110],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10110,ub_10110])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_tot_AX_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_AllCodes')

figure(1020)
lb_1020=0;
ub_1020=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_woWare_e+j_AX_woWare_d)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_1020,ub_1020],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_1020,ub_1020])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','total','Location','South')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('ambcond_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_NEO2')

figure(10200)
lb_10200=0;
ub_10200=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_woWare_e_tot+j_AX_woWare_d_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_10200,ub_10200],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_10200,ub_10200])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','total','Location','SouthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('ambcond_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_NEO2')

figure(1021)
lb_1021=0;
ub_1021=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_e+j_AX_d)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_1021,ub_1021],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_1021,ub_1021])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','total')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('ambcond_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_NEO2')

figure(10210)
lb_10210=0;
ub_10210=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_e_tot+j_AX_d_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_10210,ub_10210],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_10210,ub_10210])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','total','Location','SouthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('ambcond_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_NEO2')

%% 2 species NA results (3D VMEC equilibrium 32138 + Profiles 32169 w/ hel. core) 

ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

% surface area for computation of total fluxes
surf_area=surf_area_NEO2{ind_data_NEO2};

% Get NEO-2 data for plots
num_radial_pts=numel(Gamma_NA_spec_NEO2{ind_data_NEO2});
Gamma_NA_e=zeros(num_radial_pts,1);
Gamma_NA_d=zeros(num_radial_pts,1);
Gamma_NA_Ware_e=zeros(num_radial_pts,1);
Gamma_NA_Ware_d=zeros(num_radial_pts,1);
Qflux_NA_e=zeros(num_radial_pts,1);
Qflux_NA_d=zeros(num_radial_pts,1);
Qflux_NA_Ware_e=zeros(num_radial_pts,1);
Qflux_NA_Ware_d=zeros(num_radial_pts,1);
ne=zeros(num_radial_pts,1);
nd=zeros(num_radial_pts,1);
Te=zeros(num_radial_pts,1);
Td=zeros(num_radial_pts,1);
j_NA_e=zeros(num_radial_pts,1);
j_NA_d=zeros(num_radial_pts,1);
j_NA_Ware_e=zeros(num_radial_pts,1);
j_NA_Ware_d=zeros(num_radial_pts,1);
for k=1:num_radial_pts
    Gamma_NA_e(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_NA_d(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_NA_Ware_e(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_NA_Ware_d(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_NA_e(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_NA_d(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_NA_Ware_e(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_NA_Ware_d(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    ne(k)=n_spec_NEO2{ind_data_NEO2}{k}(1);
    nd(k)=n_spec_NEO2{ind_data_NEO2}{k}(2);
    Te(k)=T_spec_NEO2{ind_data_NEO2}{k}(1);
    Td(k)=T_spec_NEO2{ind_data_NEO2}{k}(2);
    j_NA_e(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_NA_d(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_NA_Ware_e(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_NA_Ware_d(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
end

Qflux_NA_e_tot=Qflux_NA_e.*surf_area;
Qflux_NA_d_tot=Qflux_NA_d.*surf_area;

j_NA_e_tot=j_NA_e.*surf_area;
j_NA_d_tot=j_NA_d.*surf_area;

Gamma_NA_woWare_e=Gamma_NA_e-Gamma_NA_Ware_e;
Gamma_NA_woWare_d=Gamma_NA_d-Gamma_NA_Ware_d;

Qflux_NA_woWare_e=Qflux_NA_e-Qflux_NA_Ware_e;
Qflux_NA_woWare_d=Qflux_NA_d-Qflux_NA_Ware_d;
Qflux_NA_woWare_e_tot=Qflux_NA_woWare_e.*surf_area;
Qflux_NA_woWare_d_tot=Qflux_NA_woWare_d.*surf_area;

j_NA_woWare_e=j_NA_e-j_NA_Ware_e;
j_NA_woWare_d=j_NA_d-j_NA_Ware_d;
j_NA_woWare_e_tot=j_NA_woWare_e.*surf_area;
j_NA_woWare_d_tot=j_NA_woWare_d.*surf_area;

figure(2010)
lb_2010=0;
ub_2010=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_woWare_d./nd),'r:>')  % conversion SI units
%plot([lb_2010,ub_2010],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_2010,ub_2010])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_NA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_AllCodes')

figure(20100)
lb_20100=0;
ub_20100=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_NA_woWare_e_tot,'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_NA_woWare_d_tot,'r:>')  % conversion SI units
%plot([lb_20100,ub_20100],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20100,ub_20100])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_tot_NA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_AllCodes')

figure(2011)
lb_2011=0;
ub_2011=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_d./nd),'r:>')  % conversion SI units
%plot([lb_2011,ub_2011],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_2011,ub_2011])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_NA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_AllCodes')

figure(20110)
lb_20110=0;
ub_20110=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_NA_e_tot,'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*Qflux_NA_d_tot,'r:>')  % conversion SI units
%plot([lb_20110,ub_20110],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20110,ub_20110])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('Qflux_tot_NA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_AllCodes')

figure(1020)
lb_1020=0;
ub_1020=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_NA_woWare_e+j_NA_woWare_d)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_1020,ub_1020],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_1020,ub_1020])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('jrNA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_NEO2')

figure(20200)
lb_20200=0;
ub_20200=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_NA_woWare_e_tot+j_NA_woWare_d_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_20200,ub_20200],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20200,ub_20200])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('jrNA_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_woWare_Spec2_NEO2')

figure(2021)
lb_2021=0;
ub_2021=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_NA_e+j_NA_d)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_2021,ub_2021],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_2021,ub_2021])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('jrNA_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_NEO2')

figure(20210)
lb_20210=0;
ub_20210=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_NA_e_tot+j_NA_d_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_20210,ub_20210],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20210,ub_20210])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('jrNA_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_wWare_Spec2_NEO2')

figure(3000)
lb_3000=0;
ub_3000=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),0.1*TphiNA_tot_NEO2{ind_data_NEO2},...
     'color',[0.1 0.5 0.5])
hold on
plot([lb_3000,ub_3000],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_3000,ub_3000])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$T^{\rm NA}_{\varphi} \ [{\rm Nm/m^3}]$','interpreter','latex')
legend('total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('TphiNA_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_Spec2_NEO2')

figure(3000)
lb_3000=0;
ub_3000=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-7*TphiNA_int_tot_NEO2{ind_data_NEO2},...
     'color',[0.1 0.5 0.5])
hold on
plot([lb_3000,ub_3000],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_3000,ub_3000])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$T^{\rm NA, int}_{\varphi} \ [{\rm Nm}]$','interpreter','latex')
legend('total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.0s'})
export_plot('TphiNAint_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.0s_Spec2_NEO2')

