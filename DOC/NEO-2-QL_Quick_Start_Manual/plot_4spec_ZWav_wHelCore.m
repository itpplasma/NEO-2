%% define constants
e=4.8032e-10; % elementary charge

%% load NEO-2 output
% 2 species benchmarking configuration (ExB only)
data_NEO2_nspec4_lag12=h52struct('final_neo2_multispecies_out.h5');
fname_NEO2_nspec4_lag12=fieldnames(data_NEO2_nspec4_lag12);

% allocate storage array
num_data_NEO2=1; % total number of HDFF5 data files to be processed
data_NEO2=cell(num_data_NEO2,2);
% 2 species benchmarking configuration (ExB only)
data_NEO2{1,1}=data_NEO2_nspec4_lag12;
data_NEO2{1,2}=fname_NEO2_nspec4_lag12;

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
av_nabla_stor_NEO2=cell(num_data_NEO2,1);
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
    av_nabla_stor_NEO2_tmp=zeros(numel(data_fname),1);
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
        av_nabla_stor_NEO2_tmp(es_ind)=data_struct.(data_fname{es_ind}).av_nabla_stor;
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
    av_nabla_stor_NEO2{file_ind}=av_nabla_stor_NEO2_tmp(1:data_ctr);
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
    av_nabla_stor_NEO2_tmp=av_nabla_stor_NEO2{file_ind};
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
    surf_area_NEO2{file_ind}=(4*(pi^2)*av_nabla_stor_NEO2_tmp./avb2_NEO2_tmp).*...
        abs(boozer_psi_pr_NEO2_tmp.*(aiota_NEO2_tmp.*bcovar_tht_NEO2_tmp+bcovar_phi_NEO2_tmp));
end

%return
%% Compute radial and inductive electric field for Hakan Smith

% 2 species results (3D VMEC equilibrium 31021 + Profiles 32169 wo hel. core)
ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

E_rhotor=2*Er_NEO2{ind_data_NEO2}.*sqrt(boozer_s_NEO2{ind_data_NEO2})./av_nabla_stor_NEO2{ind_data_NEO2}; % [statV]
E_rhotor=E_rhotor*3e2*1e-3; % [kV]

E_rhotor_for_Hakan=[sqrt(boozer_s_NEO2{ind_data_NEO2}),E_rhotor];

avEparB_ov_SqrtAvB2=avEparB_ov_avb2_NEO2{ind_data_NEO2}.*sqrt(avb2_NEO2{ind_data_NEO2});
avEparB_ov_SqrtAvB2=avEparB_ov_SqrtAvB2*3e4*1e3; % [mV/m]

EparB_for_Hakan=[sqrt(boozer_s_NEO2{ind_data_NEO2}),avEparB_ov_SqrtAvB2]

figure(2)
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),avEparB_ov_SqrtAvB2,'-o')
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\langle E_{\parallel} B \rangle \ / \ \langle B^2 \rangle^{1/2}$ [mV/m]',...
    'interpreter','latex')
set(gca,'XTick',[0:0.2:1])
export_plot('avEparB_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_Spec4_ZWav')

%return
%% 4 species results (3D VMEC equilibrium 31021 + Profiles 32169 w/o hel. core) 

ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

% surface area for computation of total fluxes
surf_area=surf_area_NEO2{ind_data_NEO2};

% Get NEO-2 data for plots
num_radial_pts=numel(Gamma_AX_spec_NEO2{ind_data_NEO2});
Gamma_AX_e=zeros(num_radial_pts,1);
Gamma_AX_d=zeros(num_radial_pts,1);
Gamma_AX_B=zeros(num_radial_pts,1);
Gamma_AX_W=zeros(num_radial_pts,1);
Gamma_AX_Ware_e=zeros(num_radial_pts,1);
Gamma_AX_Ware_d=zeros(num_radial_pts,1);
Gamma_AX_Ware_B=zeros(num_radial_pts,1);
Gamma_AX_Ware_W=zeros(num_radial_pts,1);
Qflux_AX_e=zeros(num_radial_pts,1);
Qflux_AX_d=zeros(num_radial_pts,1);
Qflux_AX_B=zeros(num_radial_pts,1);
Qflux_AX_W=zeros(num_radial_pts,1);
Qflux_AX_Ware_e=zeros(num_radial_pts,1);
Qflux_AX_Ware_d=zeros(num_radial_pts,1);
Qflux_AX_Ware_B=zeros(num_radial_pts,1);
Qflux_AX_Ware_W=zeros(num_radial_pts,1);
ne=zeros(num_radial_pts,1);
nd=zeros(num_radial_pts,1);
nB=zeros(num_radial_pts,1);
nW=zeros(num_radial_pts,1);
Te=zeros(num_radial_pts,1);
Td=zeros(num_radial_pts,1);
TB=zeros(num_radial_pts,1);
TW=zeros(num_radial_pts,1);
j_AX_e=zeros(num_radial_pts,1);
j_AX_d=zeros(num_radial_pts,1);
j_AX_B=zeros(num_radial_pts,1);
j_AX_W=zeros(num_radial_pts,1);
j_AX_Ware_e=zeros(num_radial_pts,1);
j_AX_Ware_d=zeros(num_radial_pts,1);
j_AX_Ware_B=zeros(num_radial_pts,1);
j_AX_Ware_W=zeros(num_radial_pts,1);
for k=1:num_radial_pts
    Gamma_AX_e(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_AX_d(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_AX_B(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(3);
    Gamma_AX_W(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(4);
    Gamma_AX_Ware_e(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_AX_Ware_d(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_AX_Ware_B(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(3);
    Gamma_AX_Ware_W(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(4);
    Qflux_AX_e(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_AX_d(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_AX_B(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(3);
    Qflux_AX_W(k)=Qflux_AX_spec_NEO2{ind_data_NEO2}{k}(4);
    Qflux_AX_Ware_e(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_AX_Ware_d(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_AX_Ware_B(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(3);
    Qflux_AX_Ware_W(k)=Qflux_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(4);
    ne(k)=n_spec_NEO2{ind_data_NEO2}{k}(1);
    nd(k)=n_spec_NEO2{ind_data_NEO2}{k}(2);
    nB(k)=n_spec_NEO2{ind_data_NEO2}{k}(3);
    nW(k)=n_spec_NEO2{ind_data_NEO2}{k}(4);
    Te(k)=T_spec_NEO2{ind_data_NEO2}{k}(1);
    Td(k)=T_spec_NEO2{ind_data_NEO2}{k}(2);
    TB(k)=T_spec_NEO2{ind_data_NEO2}{k}(3);
    TW(k)=T_spec_NEO2{ind_data_NEO2}{k}(4);
    j_AX_e(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_AX_d(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_AX_B(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(3).*zspec_NEO2{ind_data_NEO2}{k}(3)*e;
    j_AX_W(k)=Gamma_AX_spec_NEO2{ind_data_NEO2}{k}(4).*zspec_NEO2{ind_data_NEO2}{k}(4)*e;
    j_AX_Ware_e(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_AX_Ware_d(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_AX_Ware_B(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(3).*zspec_NEO2{ind_data_NEO2}{k}(3)*e;
    j_AX_Ware_W(k)=Gamma_AX_Ware_spec_NEO2{ind_data_NEO2}{k}(4).*zspec_NEO2{ind_data_NEO2}{k}(4)*e;
end
Gamma_AX_woWare_e=Gamma_AX_e-Gamma_AX_Ware_e;
Gamma_AX_woWare_d=Gamma_AX_d-Gamma_AX_Ware_d;
Gamma_AX_woWare_B=Gamma_AX_B-Gamma_AX_Ware_B;
Gamma_AX_woWare_W=Gamma_AX_W-Gamma_AX_Ware_W;
Qflux_AX_woWare_e=Qflux_AX_e-Qflux_AX_Ware_e;
Qflux_AX_woWare_d=Qflux_AX_d-Qflux_AX_Ware_d;
Qflux_AX_woWare_B=Qflux_AX_B-Qflux_AX_Ware_B;
Qflux_AX_woWare_W=Qflux_AX_W-Qflux_AX_Ware_W;
j_AX_woWare_e=j_AX_e-j_AX_Ware_e;
j_AX_woWare_d=j_AX_d-j_AX_Ware_d;
j_AX_woWare_B=j_AX_B-j_AX_Ware_B;
j_AX_woWare_W=j_AX_W-j_AX_Ware_W;

j_AX_e_tot=j_AX_e.*surf_area;
j_AX_d_tot=j_AX_d.*surf_area;
j_AX_B_tot=j_AX_B.*surf_area;
j_AX_W_tot=j_AX_W.*surf_area;
j_AX_woWare_e_tot=j_AX_woWare_e.*surf_area;
j_AX_woWare_d_tot=j_AX_woWare_d.*surf_area;
j_AX_woWare_B_tot=j_AX_woWare_B.*surf_area;
j_AX_woWare_W_tot=j_AX_woWare_W.*surf_area;

Qflux_AX_e_tot=Qflux_AX_e.*surf_area;
Qflux_AX_d_tot=Qflux_AX_d.*surf_area;
Qflux_AX_B_tot=Qflux_AX_B.*surf_area;
Qflux_AX_W_tot=Qflux_AX_W.*surf_area;
Qflux_AX_woWare_e_tot=Qflux_AX_woWare_e.*surf_area;
Qflux_AX_woWare_d_tot=Qflux_AX_woWare_d.*surf_area;
Qflux_AX_woWare_B_tot=Qflux_AX_woWare_B.*surf_area;
Qflux_AX_woWare_W_tot=Qflux_AX_woWare_W.*surf_area;

% Get NEO data for plots
% Note G_a / n_a in units m/s and Q_a / n_a in units 10^-19 W m !!!
%NEOdata = load('AM_NEO_results_4sp_GandQ'); % by C. Angioni (14.06.17)
%NEOdata(:,6:end) = NEOdata(:,6:end)*1e-19;
% possible conversion factor for NEO fluxes - seems to be small
%r_minor_LCFS=60.147232; % [cm] 
%conv_fac_reff=r_minor_LCFS*av_nabla_stor_NEO2{ind_data_NEO2}./(2*sqrt(boozer_s_NEO2{ind_data_NEO2}));
%conv_fac_reff_interp=interp1(sqrt(boozer_s_NEO2{ind_data_NEO2}),conv_fac_reff,NEOdata(:,1),'linear','extrap');

% Get NCLASS data for plots
% Note G_a / n_a in units m/s and Q_a / n_a in units 10^-19 W m !!!
%NCLASSdata = load('AM_NCL_results_4sp_GandQ'); % by C. Angioni (14.06.17)
%NCLASSdata(:,6:end) = NCLASSdata(:,6:end)*1e-19;

% Get SFINCS data for plots
%SFINCSdata = load('SFINCS_results_2sp.dat'); % by H. Smith (04.07.17)
%SFINCSdata = load('SFINCS_results_2sp_CoulLog.dat'); % by H. Smith (04.07.17)

figure(1000)
lb_1000=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
rb_1000=0.9;
tb_1000=0.03;
bb_1000=0.0;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_d./nd),'r:>')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,5),'b--')
% plot(NEOdata(:,1),NEOdata(:,2),'r--')
% plot(NCLASSdata(:,1),NCLASSdata(:,5),'b-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,2),'r-.')
% plot(SFINCSdata(:,1),SFINCSdata(:,2),'b:') 
% plot(SFINCSdata(:,1),SFINCSdata(:,3),'r:') 
plot([lb_1000,rb_1000],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1000,rb_1000])
% ylim([bb_1000,tb_1000])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
% legend('NEO: e','NEO: d','NCLASS: e','NCLASS: d','NEO-2: e','NEO-2: d','Location','EastOutside')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Gamma_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes')

figure(10001)
lb_10001=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
rb_10001=0.9;
tb_10001=0.5;
bb_10001=-2;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_B./nB),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_W./nW),'k:sq')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,3),'m--')
% plot(NEOdata(:,1),NEOdata(:,4),'k--')
% plot(NCLASSdata(:,1),NCLASSdata(:,3),'m-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,4),'k-.')
plot([lb_10001,rb_10001],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10001,rb_10001])
ylim([bb_10001,tb_10001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
% legend('NEO: B','NEO: W','NCLASS: B','NCLASS: W','NEO-2: B','NEO-2: W','Location','EastOutside')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Gamma_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes_zoom')

%return

figure(1001)
lb_1001=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
rb_1001=1.2;
tb_1001=0.03;
bb_1001=-0.06;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_e./ne),'b:o') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_d./nd),'r:+')  % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_e./ne),'b--sq') % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_d./nd),'r-->')  % conversion SI units
plot([lb_1001,rb_1001],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1001,rb_1001])
%ylim([bb_1001,tb_1001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
lh=legend('e: w/ $E_\parallel$','d: w/ $E_\parallel$',...
    'e: w/o $E_\parallel$','d: w/o $E_\parallel$',...
    'Location','SouthEast');
set(lh,'interpreter','latex')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Gamma_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15_woWare_Spec4_ZWav_NEO2')

figure(10010)
lb_1001=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
rb_1001=1.2;
tb_1001=0.03;
bb_1001=-0.06;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_B./nB),'m:x')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_W./nW),'k:*')  % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_B./nB),'m--d')  % conversion SI units
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-2*(Gamma_AX_woWare_W./nW),'k--h')  % conversion SI units
plot([lb_1001,rb_1001],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_1001,rb_1001])
%ylim([bb_1001,tb_1001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$\Gamma^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm m s^{-1}}]$','interpreter','latex')
lh=legend('B: w/ $E_\parallel$','W: w/ $E_\parallel$',...
    'B: w/o $E_\parallel$','W: w/o $E_\parallel$','Location','SouthEast');
set(lh,'interpreter','latex')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Gamma_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2_zoom')

%return

figure(1010)
lb_1010=0;
ub_1010=0.9;
%tb_1010=1e-15;
%bb_1010=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-9*(Qflux_AX_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-9*(Qflux_AX_woWare_d./nd),'r:>')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,9),'b--')
% plot(NEOdata(:,1),NEOdata(:,6),'r--')
% plot(NCLASSdata(:,1),NCLASSdata(:,9),'b-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,6),'r-.')
%plot([lb_1010,ub_1010],[0,0],'color',[0.8 0.8 0.8]);
hold off
%set(gca,'YScale','log')
xlim([lb_1010,ub_1010])
%ylim([bb_1010,tb_1010])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [{\rm W m}]$','interpreter','latex')
% legend('NEO: e','NEO: d','NCLASS: e','NCLASS: d','NEO-2: e','NEO-2: d',...
%     'Location','EastOutside')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes')

figure(1011)
lb_1011=0;
ub_1011=0.9;
%tb_1011=1e-15;
%bb_1011=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_AX_woWare_e_tot),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_AX_woWare_d_tot),'r:>')  % conversion SI units
%plot([lb_1011,ub_1011],[0,0],'color',[0.8 0.8 0.8]);
hold off
%set(gca,'YScale','log')
xlim([lb_1011,ub_1011])
%ylim([bb_1011,tb_1011])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes')

figure(1012)
lb_1012=0;
ub_1012=0.9;
%tb_1012=1e-15;
%bb_1012=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_AX_d./nd),'r:>')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,9),'b--')
% plot(NEOdata(:,1),NEOdata(:,6),'r--')
% plot(NCLASSdata(:,1),NCLASSdata(:,9),'b-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,6),'r-.')
%plot([lb_1012,ub_1012],[0,0],'color',[0.8 0.8 0.8]);
hold off
set(gca,'YScale','linear')
xlim([lb_1012,ub_1012])
%ylim([bb_1012,tb_1012])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
% legend('NEO: e','NEO: d','NCLASS: e','NCLASS: d','NEO-2: e','NEO-2: d',...
%     'Location','EastOutside')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes')

figure(1013)
lb_1013=0;
ub_1013=0.9;
%tb_1013=1e-15;
%bb_1013=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_AX_e_tot),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_AX_d_tot),'r:>')  % conversion SI units
%plot([lb_1013,ub_1013],[0,0],'color',[0.8 0.8 0.8]);
hold off
set(gca,'YScale','linear')
xlim([lb_1013,ub_1013])
%ylim([bb_1013,tb_1013])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes')

figure(10100)
lb_10100=0;
ub_10100=0.9;
%tb_10100=1e-15;
%bb_10100=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_AX_woWare_B./nB),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_AX_woWare_W./nW),'k:sq')  % conversion SI units
% plot(NEOdata(:,1),1e15*NEOdata(:,7),'m--')
% plot(NEOdata(:,1),1e15*NEOdata(:,8),'k--')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,7),'m-.')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,8),'k-.')
plot([lb_10100,ub_10100],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10100,ub_10100])
%ylim([bb_10100,tb_10100])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [10^{-15} \ {\rm W m}]$','interpreter','latex')
%legend('NEO: B','NEO: W','NCLASS: B','NCLASS: W','NEO-2: B','NEO-2: W','Location','SouthEast')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes_zoom')

figure(10101)
lb_10101=0;
ub_10101=0.9;
%tb_10101=1e-15;
%bb_10101=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_AX_woWare_B_tot),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_AX_woWare_W_tot),'k:sq')  % conversion SI units
plot([lb_10101,ub_10101],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10101,ub_10101])
%ylim([bb_10101,tb_10101])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm kW}]$','interpreter','latex')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes_zoom')

figure(10102)
lb_10102=0;
ub_10102=0.9;
%tb_10102=1e-15;
%bb_10102=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_AX_B./nB),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_AX_W./nW),'k:sq')  % conversion SI units
% plot(NEOdata(:,1),1e15*NEOdata(:,7),'m--')
% plot(NEOdata(:,1),1e15*NEOdata(:,8),'k--')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,7),'m-.')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,8),'k-.')
plot([lb_10102,ub_10102],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10102,ub_10102])
%ylim([bb_10102,tb_10102])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha n_\alpha^{-1} \ [10^{-15} \ {\rm W m}]$','interpreter','latex')
%legend('NEO: B','NEO: W','NCLASS: B','NCLASS: W','NEO-2: B','NEO-2: W','Location','SouthEast')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes_zoom')

figure(10103)
lb_10103=0;
ub_10103=0.9;
%tb_10103=1e-15;
%bb_10103=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_AX_B_tot),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_AX_W_tot),'k:sq')  % conversion SI units
plot([lb_10103,ub_10103],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_10103,ub_10103])
%ylim([bb_10103,tb_10103])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm AX}_\alpha \ [{\rm kW}]$','interpreter','latex')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_AX_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes_zoom')

%return

figure(1020)
lb_1020=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_1020=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_B/3e5,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_W/3e5,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_woWare_e+j_AX_woWare_d+j_AX_woWare_B+j_AX_woWare_W)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_1020,ub_1020],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_1020,ub_1020])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','B','W','total','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('ambcond_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2')

figure(10200)
lb_10200=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_10200=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_B_tot/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_woWare_W_tot/3e9,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_woWare_e_tot+j_AX_woWare_d_tot+j_AX_woWare_B_tot+j_AX_woWare_W_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_10200,ub_10200],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_10200,ub_10200])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','total','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('ambcond_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2')

figure(1021)
lb_1021=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_1021=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_e/3e5,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_d/3e5,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_B/3e5,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_W/3e5,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_e+j_AX_d+j_AX_B+j_AX_W)/3e5,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_1021,ub_1021],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_1021,ub_1021])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A m^{-2}}]$','interpreter','latex')
legend('e','d','B','W','total')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('ambcond_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2')

figure(10210)
lb_10210=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_10210=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_B_tot/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_AX_W_tot/3e9,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_AX_e_tot+j_AX_d_tot+j_AX_B_tot+j_AX_W_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_10210,ub_10210],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_10210,ub_10210])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm AX}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','total')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('ambcond_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2')

%return

% save NEO-2 data for Hakan Smith and Clemente Angioni
% -> NEO-2 data
empty_col=zeros(size(boozer_s_NEO2{ind_data_NEO2}));
Mout_NEO2=[...
    sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    1e-2*(Gamma_AX_woWare_e./ne),...
    1e-2*(Gamma_AX_woWare_d./nd),...
    1e-2*(Gamma_AX_woWare_B./nB),...
    1e-2*(Gamma_AX_woWare_W./nW),...
    1e-9*(Qflux_AX_woWare_e./ne),...
    1e-9*(Qflux_AX_woWare_d./nd),...
    1e-9*(Qflux_AX_woWare_B./nB),...
    1e-9*(Qflux_AX_woWare_W./nW)];
[~,sortInd]=sort(Mout_NEO2(:,1));
Mout_NEO2=Mout_NEO2(sortInd,:);
% -> file format
formatSep='   ';
col1HdrLabel=['rho_tor',repmat(' ',1,6)];
col2HdrLabel=['Ge/ne [m/s]',repmat(' ',1,2)];
col3HdrLabel=['G1/n1 [m/s]',repmat(' ',1,2)];
col4HdrLabel=['G2/n2 [m/s]',repmat(' ',1,2)];
col5HdrLabel=['G3/n3 [m/s]',repmat(' ',1,2)];
col6HdrLabel=['Qe/ne [W m]',repmat(' ',1,2)];
col7HdrLabel=['Q1/n1 [W m]',repmat(' ',1,2)];
col8HdrLabel=['Q2/n2 [W m]',repmat(' ',1,2)];
col9HdrLabel=['Q3/n3 [W m]',repmat(' ',1,2)];
formatSpecHdr=[repmat([formatSep,'%s'],1,9),'\n'];
formatSpecHdr=['%%',formatSpecHdr(2:end)];
formatSepData='  ';
formatSpecData=[repmat([formatSepData,'% 8.7e'],1,9),'\n'];
% -> write file
fileID=fopen('NEO2_results_4sp_GandQ_PROF_AUG32169_at_t=4.15s.dat','w');
fprintf(fileID,formatSpecHdr,col1HdrLabel,...
    col2HdrLabel,col3HdrLabel,col4HdrLabel,col5HdrLabel,...
    col6HdrLabel,col7HdrLabel,col8HdrLabel,col9HdrLabel);
for line_ptr=1:size(Mout_NEO2,1)
    fprintf(fileID,formatSpecData,Mout_NEO2(line_ptr,:));
end
fclose(fileID);
% % display data
% str_sep=repmat('  ',numel(Mout_NEO2(:,1)),1);
% Mout2str_NEO2=[...
%     num2str(Mout_NEO2(:,1),'%5.4f'),str_sep,...
%     num2str(Mout_NEO2(:,2),'%5.4e'),str_sep,...
%     num2str(Mout_NEO2(:,3),'%5.4e')];
% disp('rho_tor Flux_e      Flux_D')
% disp(Mout2str_NEO2)
 
%return

%% %% 4 species NA results (3D VMEC equilibrium 31021 + Profiles 32169 w/o hel. core) 

ind_data_NEO2=1; % select NEO-2 output data file (HDF5)

% surface area for computation of total fluxes
surf_area=surf_area_NEO2{ind_data_NEO2};

% Get NEO-2 data for plots
num_radial_pts=numel(Gamma_NA_spec_NEO2{ind_data_NEO2});
Gamma_NA_e=zeros(num_radial_pts,1);
Gamma_NA_d=zeros(num_radial_pts,1);
Gamma_NA_B=zeros(num_radial_pts,1);
Gamma_NA_W=zeros(num_radial_pts,1);
Gamma_NA_Ware_e=zeros(num_radial_pts,1);
Gamma_NA_Ware_d=zeros(num_radial_pts,1);
Gamma_NA_Ware_B=zeros(num_radial_pts,1);
Gamma_NA_Ware_W=zeros(num_radial_pts,1);
Qflux_NA_e=zeros(num_radial_pts,1);
Qflux_NA_d=zeros(num_radial_pts,1);
Qflux_NA_B=zeros(num_radial_pts,1);
Qflux_NA_W=zeros(num_radial_pts,1);
Qflux_NA_Ware_e=zeros(num_radial_pts,1);
Qflux_NA_Ware_d=zeros(num_radial_pts,1);
Qflux_NA_Ware_B=zeros(num_radial_pts,1);
Qflux_NA_Ware_W=zeros(num_radial_pts,1);
ne=zeros(num_radial_pts,1);
nd=zeros(num_radial_pts,1);
nB=zeros(num_radial_pts,1);
nW=zeros(num_radial_pts,1);
Te=zeros(num_radial_pts,1);
Td=zeros(num_radial_pts,1);
TB=zeros(num_radial_pts,1);
TW=zeros(num_radial_pts,1);
j_NA_e=zeros(num_radial_pts,1);
j_NA_d=zeros(num_radial_pts,1);
j_NA_B=zeros(num_radial_pts,1);
j_NA_W=zeros(num_radial_pts,1);
j_NA_Ware_e=zeros(num_radial_pts,1);
j_NA_Ware_d=zeros(num_radial_pts,1);
j_NA_Ware_B=zeros(num_radial_pts,1);
j_NA_Ware_W=zeros(num_radial_pts,1);
for k=1:num_radial_pts
    Gamma_NA_e(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_NA_d(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_NA_B(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(3);
    Gamma_NA_W(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(4);
    Gamma_NA_Ware_e(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Gamma_NA_Ware_d(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Gamma_NA_Ware_B(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(3);
    Gamma_NA_Ware_W(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(4);
    Qflux_NA_e(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_NA_d(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_NA_B(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(3);
    Qflux_NA_W(k)=Qflux_NA_spec_NEO2{ind_data_NEO2}{k}(4);
    Qflux_NA_Ware_e(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1);
    Qflux_NA_Ware_d(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2);
    Qflux_NA_Ware_B(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(3);
    Qflux_NA_Ware_W(k)=Qflux_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(4);
    ne(k)=n_spec_NEO2{ind_data_NEO2}{k}(1);
    nd(k)=n_spec_NEO2{ind_data_NEO2}{k}(2);
    nB(k)=n_spec_NEO2{ind_data_NEO2}{k}(3);
    nW(k)=n_spec_NEO2{ind_data_NEO2}{k}(4);
    Te(k)=T_spec_NEO2{ind_data_NEO2}{k}(1);
    Td(k)=T_spec_NEO2{ind_data_NEO2}{k}(2);
    TB(k)=T_spec_NEO2{ind_data_NEO2}{k}(3);
    TW(k)=T_spec_NEO2{ind_data_NEO2}{k}(4);
    j_NA_e(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_NA_d(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_NA_B(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(3).*zspec_NEO2{ind_data_NEO2}{k}(3)*e;
    j_NA_W(k)=Gamma_NA_spec_NEO2{ind_data_NEO2}{k}(4).*zspec_NEO2{ind_data_NEO2}{k}(4)*e;
    j_NA_Ware_e(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(1).*zspec_NEO2{ind_data_NEO2}{k}(1)*e;
    j_NA_Ware_d(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(2).*zspec_NEO2{ind_data_NEO2}{k}(2)*e;
    j_NA_Ware_B(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(3).*zspec_NEO2{ind_data_NEO2}{k}(3)*e;
    j_NA_Ware_W(k)=Gamma_NA_Ware_spec_NEO2{ind_data_NEO2}{k}(4).*zspec_NEO2{ind_data_NEO2}{k}(4)*e;
end
Gamma_NA_woWare_e=Gamma_NA_e-Gamma_NA_Ware_e;
Gamma_NA_woWare_d=Gamma_NA_d-Gamma_NA_Ware_d;
Gamma_NA_woWare_B=Gamma_NA_B-Gamma_NA_Ware_B;
Gamma_NA_woWare_W=Gamma_NA_W-Gamma_NA_Ware_W;
Qflux_NA_woWare_e=Qflux_NA_e-Qflux_NA_Ware_e;
Qflux_NA_woWare_d=Qflux_NA_d-Qflux_NA_Ware_d;
Qflux_NA_woWare_B=Qflux_NA_B-Qflux_NA_Ware_B;
Qflux_NA_woWare_W=Qflux_NA_W-Qflux_NA_Ware_W;
j_NA_woWare_e=j_NA_e-j_NA_Ware_e;
j_NA_woWare_d=j_NA_d-j_NA_Ware_d;
j_NA_woWare_B=j_NA_B-j_NA_Ware_B;
j_NA_woWare_W=j_NA_W-j_NA_Ware_W;

j_NA_e_tot=j_NA_e.*surf_area;
j_NA_d_tot=j_NA_d.*surf_area;
j_NA_B_tot=j_NA_B.*surf_area;
j_NA_W_tot=j_NA_W.*surf_area;
j_NA_woWare_e_tot=j_NA_woWare_e.*surf_area;
j_NA_woWare_d_tot=j_NA_woWare_d.*surf_area;
j_NA_woWare_B_tot=j_NA_woWare_B.*surf_area;
j_NA_woWare_W_tot=j_NA_woWare_W.*surf_area;

Qflux_NA_e_tot=Qflux_NA_e.*surf_area;
Qflux_NA_d_tot=Qflux_NA_d.*surf_area;
Qflux_NA_B_tot=Qflux_NA_B.*surf_area;
Qflux_NA_W_tot=Qflux_NA_W.*surf_area;
Qflux_NA_woWare_e_tot=Qflux_NA_woWare_e.*surf_area;
Qflux_NA_woWare_d_tot=Qflux_NA_woWare_d.*surf_area;
Qflux_NA_woWare_B_tot=Qflux_NA_woWare_B.*surf_area;
Qflux_NA_woWare_W_tot=Qflux_NA_woWare_W.*surf_area;

figure(2010)
lb_2010=0;
ub_2010=0.9;
%tb_2010=1e-15;
%bb_2010=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-9*(Qflux_NA_woWare_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-9*(Qflux_NA_woWare_d./nd),'r:>')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,9),'b--')
% plot(NEOdata(:,1),NEOdata(:,6),'r--')
% plot(NCLASSdata(:,1),NCLASSdata(:,9),'b-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,6),'r-.')
%plot([lb_2010,ub_2010],[0,0],'color',[0.8 0.8 0.8]);
hold off
%set(gca,'YScale','log')
xlim([lb_2010,ub_2010])
%ylim([bb_2010,tb_2010])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [{\rm W m}]$','interpreter','latex')
% legend('NEO: e','NEO: d','NCLASS: e','NCLASS: d','NEO-2: e','NEO-2: d',...
%     'Location','EastOutside')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes')

figure(2011)
lb_2011=0;
ub_2011=0.9;
%tb_2011=1e-15;
%bb_2011=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_NA_woWare_e_tot),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_NA_woWare_d_tot),'r:>')  % conversion SI units
%plot([lb_2011,ub_2011],[0,0],'color',[0.8 0.8 0.8]);
hold off
%set(gca,'YScale','log')
xlim([lb_2011,ub_2011])
%ylim([bb_2011,tb_2011])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes')

figure(2012)
lb_2012=0;
ub_2012=0.9;
%tb_2012=1e-15;
%bb_2012=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_e./ne),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e16*1e-9*(Qflux_NA_d./nd),'r:>')  % conversion SI units
% plot(NEOdata(:,1),NEOdata(:,9),'b--')
% plot(NEOdata(:,1),NEOdata(:,6),'r--')
% plot(NCLASSdata(:,1),NCLASSdata(:,9),'b-.')
% plot(NCLASSdata(:,1),NCLASSdata(:,6),'r-.')
%plot([lb_2012,ub_2012],[0,0],'color',[0.8 0.8 0.8]);
hold off
set(gca,'YScale','linear')
xlim([lb_2012,ub_2012])
%ylim([bb_2012,tb_2012])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [10^{-16} \ {\rm W m}]$','interpreter','latex')
% legend('NEO: e','NEO: d','NCLASS: e','NCLASS: d','NEO-2: e','NEO-2: d',...
%     'Location','EastOutside')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes')

figure(2013)
lb_2013=0;
ub_2013=0.9;
%tb_2013=1e-15;
%bb_2013=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_NA_e_tot),'b:sq') % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-6*1e-7*(Qflux_NA_d_tot),'r:>')  % conversion SI units
%plot([lb_2013,ub_2013],[0,0],'color',[0.8 0.8 0.8]);
hold off
set(gca,'YScale','linear')
xlim([lb_2013,ub_2013])
%ylim([bb_2013,tb_2013])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm MW}]$','interpreter','latex')
legend('NEO-2: e','NEO-2: d','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes')

figure(20100)
lb_20100=0;
ub_20100=0.9;
%tb_20100=1e-15;
%bb_20100=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_NA_woWare_B./nB),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_NA_woWare_W./nW),'k:sq')  % conversion SI units
% plot(NEOdata(:,1),1e15*NEOdata(:,7),'m--')
% plot(NEOdata(:,1),1e15*NEOdata(:,8),'k--')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,7),'m-.')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,8),'k-.')
plot([lb_20100,ub_20100],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20100,ub_20100])
%ylim([bb_20100,tb_20100])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [10^{-15} \ {\rm W m}]$','interpreter','latex')
%legend('NEO: B','NEO: W','NCLASS: B','NCLASS: W','NEO-2: B','NEO-2: W','Location','SouthEast')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes_zoom')

figure(20101)
lb_20101=0;
ub_20101=0.9;
%tb_20101=1e-15;
%bb_20101=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_NA_woWare_B_tot),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_NA_woWare_W_tot),'k:sq')  % conversion SI units
plot([lb_20101,ub_20101],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20101,ub_20101])
%ylim([bb_20101,tb_20101])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm kW}]$','interpreter','latex')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_AllCodes_zoom')

figure(20102)
lb_20102=0;
ub_20102=0.9;
%tb_20102=1e-15;
%bb_20102=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_NA_B./nB),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e15*1e-9*(Qflux_NA_W./nW),'k:sq')  % conversion SI units
% plot(NEOdata(:,1),1e15*NEOdata(:,7),'m--')
% plot(NEOdata(:,1),1e15*NEOdata(:,8),'k--')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,7),'m-.')
% plot(NCLASSdata(:,1),1e15*NCLASSdata(:,8),'k-.')
plot([lb_20102,ub_20102],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20102,ub_20102])
%ylim([bb_20102,tb_20102])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha n_\alpha^{-1} \ [10^{-15} \ {\rm W m}]$','interpreter','latex')
%legend('NEO: B','NEO: W','NCLASS: B','NCLASS: W','NEO-2: B','NEO-2: W','Location','SouthEast')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes_zoom')

figure(20103)
lb_20103=0;
ub_20103=0.9;
%tb_20103=1e-15;
%bb_20103=1e-18;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_NA_B_tot),'m:d')  % conversion SI units
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-3*1e-7*(Qflux_NA_W_tot),'k:sq')  % conversion SI units
plot([lb_20103,ub_20103],[0,0],'color',[0.8 0.8 0.8]);
hold off
xlim([lb_20103,ub_20103])
%ylim([bb_20103,tb_20103])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$Q^{\rm NA}_\alpha \ [{\rm kW}]$','interpreter','latex')
legend('NEO-2: B','NEO-2: W','Location','SouthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('Qflux_NA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_AllCodes_zoom')

figure(20200)
lb_20200=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_20200=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_B_tot/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_W_tot/3e9,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_NA_woWare_e_tot+j_NA_woWare_d_tot+j_NA_woWare_B_tot+j_NA_woWare_W_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_20200,ub_20200],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20200,ub_20200])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','total','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jrNA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2')

figure(202001)
lb_202001=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_202001=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_B_tot/3e9,'m:x')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_woWare_W_tot/3e9,'k:sq')
plot([lb_202001,ub_202001],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_202001,ub_202001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('B','W','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jrNA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2_zoom')

figure(20201)
lb_20201=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_20201=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_e_tot+j_NA_woWare_e_tot)/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_d_tot+j_NA_woWare_d_tot)/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_B_tot+j_NA_woWare_B_tot)/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_W_tot+j_NA_woWare_W_tot)/3e9,'k:sq')
plot([lb_20201,ub_20201],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20201,ub_20201])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jr_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2')

figure(202011)
lb_202011=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_202011=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_B_tot+j_NA_woWare_B_tot)/3e9,'m:x')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_woWare_W_tot+j_NA_woWare_W_tot)/3e9,'k:sq')
plot([lb_202011,ub_202011],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_202011,ub_202011])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r} \ [{\rm A}]$','interpreter','latex')
legend('B','W','Location','NorthEast')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jr_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_woWare_Spec4_ZWav_NEO2_zoom')

figure(20210)
lb_20210=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_20210=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_e_tot/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_d_tot/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_B_tot/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_W_tot/3e9,'k:sq')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),...
    (j_NA_e_tot+j_NA_d_tot+j_NA_B_tot+j_NA_W_tot)/3e9,...
    'color',[0.1 0.5 0.5],'marker','>');%,'linewidth',2)
plot([lb_20210,ub_20210],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20210,ub_20210])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','total')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jrNA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2')

figure(202101)
lb_202101=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_202101=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_B_tot/3e9,'m:x')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),j_NA_W_tot/3e9,'k:sq')
plot([lb_202101,ub_202101],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_202101,ub_202101])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r,{\rm NA}} \ [{\rm A}]$','interpreter','latex')
legend('B','W')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jrNA_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2_zoom')

figure(20211)
lb_20211=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_20211=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_e_tot+j_NA_e_tot)/3e9,'b:o')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_d_tot+j_NA_d_tot)/3e9,'r:+')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_B_tot+j_NA_B_tot)/3e9,'m:x')
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_W_tot+j_NA_W_tot)/3e9,'k:sq')
plot([lb_20211,ub_20211],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_20211,ub_20211])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r} \ [{\rm A}]$','interpreter','latex')
legend('e','d','B','W','total')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jr_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2')

figure(202111)
lb_202111=min(sqrt(boozer_s_NEO2{ind_data_NEO2}))-2e-2;
ub_202111=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_B_tot+j_NA_B_tot)/3e9,'m:x')
hold on
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),(j_AX_W_tot+j_NA_W_tot)/3e9,'k:sq')
plot([lb_202111,ub_202111],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_202111,ub_202111])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$j^{r} \ [{\rm A}]$','interpreter','latex')
legend('B','W')
title({'VMEC: AUG #31021','PROFILE: AUG #32169 at t=4.15s'})
export_plot('jr_tot_VMEC_AUG31021_PROF_AUG32169_at_t=4.15s_wWare_Spec4_ZWav_NEO2_zoom')

figure(3000)
lb_3000=0;
ub_3000=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),0.1*TphiNA_tot_NEO2{ind_data_NEO2},...
     'color',[0.1 0.5 0.5],'marker','>')
hold on
plot([lb_3000,ub_3000],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_3000,ub_3000])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$T^{\rm NA}_{\varphi} \ [{\rm Nm/m^3}]$','interpreter','latex')
legend('total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.15s'})
export_plot('TphiNA_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.15s_Spec2_ZWav_NEO2')

figure(3001)
lb_3001=0;
ub_3001=0.9;
plot(sqrt(boozer_s_NEO2{ind_data_NEO2}),1e-7*TphiNA_int_tot_NEO2{ind_data_NEO2},...
     'color',[0.1 0.5 0.5],'marker','>')
hold on
plot([lb_3001,ub_3001],[0,0],'color',[0.8 0.8 0.8 0.5]);
hold off
xlim([lb_3001,ub_3001])
xlabel('$\rho_{\rm tor}$','interpreter','latex')
ylabel('$T^{\rm NA, int}_{\varphi} \ [{\rm Nm}]$','interpreter','latex')
legend('total','Location','NorthEast')
title({'VMEC: AUG #32138','PROFILE: AUG #32169 at t=4.15s'})
export_plot('TphiNAint_tot_VMEC_AUG32138_PROF_AUG32169_at_t=4.15s_Spec2_ZWav_NEO2')