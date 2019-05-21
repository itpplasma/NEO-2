function export_2spec_Matyas()

  [boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io, Mt_e, Mt_d] = data_process('.');
  plot_2spec_export(boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io)

end

function [boozer_s_NEO2, TphiNA_tot_NEO2, TphiNA_int_tot_NEO2, TphiNA_int_ele_NEO2, TphiNA_int_io_NEO2, Mt_e_NEO2, Mt_d_NEO2] = data_process(folder)
  %% define constants
  e=4.8032e-10; % elementary charge

  %% load NEO-2 output
  % 2 species (ExB only)
  [folder,'/final_neo2_multispecies_out.h5']
  data_NEO2_nspec2_lag7_newint_wHelCore = h52struct([folder,'/final_neo2_multispecies_out.h5']);
  fname_NEO2_nspec2_lag7_newint_wHelCore = fieldnames(data_NEO2_nspec2_lag7_newint_wHelCore);

  % allocate storage array
  num_data_NEO2 = 1; % total number of HDFF5 data files to be processed
  data_NEO2 = cell(num_data_NEO2,2);
  % 2 species (ExB only)
  data_NEO2{1,1} = data_NEO2_nspec2_lag7_newint_wHelCore;
  data_NEO2{1,2} = fname_NEO2_nspec2_lag7_newint_wHelCore;

  % extract data
  mspec_NEO2 = cell(num_data_NEO2,1);
  zspec_NEO2 = cell(num_data_NEO2,1);
  boozer_s_NEO2 = cell(num_data_NEO2,1);
  boozer_psi_pr_NEO2 = cell(num_data_NEO2,1);
  aiota_NEO2 = cell(num_data_NEO2,1);
  bcovar_tht_NEO2 = cell(num_data_NEO2,1);
  bcovar_phi_NEO2 = cell(num_data_NEO2,1);
  avb2_NEO2 = cell(num_data_NEO2,1);
  Er_NEO2 = cell(num_data_NEO2,1);
  avEparB_ov_avb2_NEO2 = cell(num_data_NEO2,1);
  avnabpsi_NEO2 = cell(num_data_NEO2,1);
  Mt_e_NEO2 = cell(num_data_NEO2,1);
  Mt_d_NEO2 = cell(num_data_NEO2,1);
  TphiNA_tot_NEO2 = cell(num_data_NEO2,1);
  TphiNA_spec_NEO2 = cell(num_data_NEO2,1);
  Gamma_AX_spec_NEO2 = cell(num_data_NEO2,1);
  Gamma_AX_Ware_spec_NEO2 = cell(num_data_NEO2,1);
  Gamma_NA_spec_NEO2 = cell(num_data_NEO2,1);
  Gamma_NA_Ware_spec_NEO2 = cell(num_data_NEO2,1);
  n_spec_NEO2 = cell(num_data_NEO2,1);
  T_spec_NEO2 = cell(num_data_NEO2,1);
  Qflux_AX_spec_NEO2 = cell(num_data_NEO2,1);
  Qflux_AX_Ware_spec_NEO2 = cell(num_data_NEO2,1);
  Qflux_NA_spec_NEO2 = cell(num_data_NEO2,1);
  Qflux_NA_Ware_spec_NEO2 = cell(num_data_NEO2,1);
  for file_ind = 1:num_data_NEO2
    % get structure
    data_struct = data_NEO2{file_ind,1};
    data_fname = data_NEO2{file_ind,2};
    % get NEO-2 data
    data_ctr = 0;
    mspec_NEO2_tmp = cell(numel(data_fname),1);
    zspec_NEO2_tmp = cell(numel(data_fname),1);
    boozer_s_NEO2_tmp = zeros(numel(data_fname),1);
    psi_pr_hat_NEO2_tmp = zeros(numel(data_fname),1);
    Bref_NEO2_tmp = zeros(numel(data_fname),1);
    aiota_NEO2_tmp = zeros(numel(data_fname),1);
    bcovar_tht_NEO2_tmp = zeros(numel(data_fname),1);
    bcovar_phi_NEO2_tmp = zeros(numel(data_fname),1);
    avbhat2_NEO2_tmp = zeros(numel(data_fname),1);
    R0_NEO2_tmp = zeros(numel(data_fname),1);
    Er_NEO2_tmp = zeros(numel(data_fname),1);
    avEparB_ov_avb2_NEO2_tmp = zeros(numel(data_fname),1);
    avnabpsi_NEO2_tmp = zeros(numel(data_fname),1);
    MteOvR_NEO2_tmp = zeros(numel(data_fname),1);
    MtdOvR_NEO2_tmp = zeros(numel(data_fname),1);
    TphiNA_tot_NEO2_tmp = zeros(numel(data_fname),1);
    TphiNA_spec_NEO2_tmp = cell(numel(data_fname),1);
    Gamma_AX_spec_NEO2_tmp = cell(numel(data_fname),1);
    Gamma_AX_Ware_spec_NEO2_tmp = cell(numel(data_fname),1);
    Gamma_NA_spec_NEO2_tmp = cell(numel(data_fname),1);
    Gamma_NA_Ware_spec_NEO2_tmp = cell(numel(data_fname),1);
    n_spec_NEO2_tmp = cell(numel(data_fname),1);
    T_spec_NEO2_tmp = cell(numel(data_fname),1);
    Qflux_AX_spec_NEO2_tmp = cell(numel(data_fname),1);
    Qflux_AX_Ware_spec_NEO2_tmp = cell(numel(data_fname),1);
    Qflux_NA_spec_NEO2_tmp = cell(numel(data_fname),1);
    Qflux_NA_Ware_spec_NEO2_tmp = cell(numel(data_fname),1);
    for es_ind = 1:numel(data_fname)
      if (~strcmp(data_fname{es_ind}(1:3),'es_'))
        continue
      end
      data_ctr = data_ctr+1;
      mspec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).m_spec;
      zspec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).z_spec;
      es_str = data_fname{es_ind}(4:end);
      boozer_s_NEO2_tmp(es_ind) = str2num(strrep(es_str,'p','.'));
      psi_pr_hat_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).psi_pr_hat;
      Bref_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).Bref;
      aiota_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).aiota;
      bcovar_tht_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).bcovar_tht;
      bcovar_phi_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).bcovar_phi;
      avbhat2_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).avbhat2;
      R0_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).R0;
      Er_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).Er;
      avEparB_ov_avb2_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).avEparB_ov_avb2;
      avnabpsi_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).avnabpsi;
      MteOvR_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).MtOvR(1);
      MtdOvR_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).MtOvR(2);
      TphiNA_tot_NEO2_tmp(es_ind) = data_struct.(data_fname{es_ind}).TphiNA_tot;
      TphiNA_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).TphiNA_spec;
      Gamma_AX_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Gamma_AX_spec;
      Gamma_AX_Ware_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Gamma_AX_Ware_spec;
      Gamma_NA_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Gamma_NA_spec;
      Gamma_NA_Ware_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Gamma_NA_Ware_spec;
      n_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).n_spec;
      T_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).T_spec;
      Qflux_AX_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Qflux_AX_spec;
      Qflux_AX_Ware_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Qflux_AX_Ware_spec;
      Qflux_NA_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Qflux_NA_spec;
      Qflux_NA_Ware_spec_NEO2_tmp{es_ind} = data_struct.(data_fname{es_ind}).Qflux_NA_Ware_spec;
    end
    % store data for plotting
    mspec_NEO2{file_ind} = mspec_NEO2_tmp(1:data_ctr);
    zspec_NEO2{file_ind} = zspec_NEO2_tmp(1:data_ctr);
    boozer_s_NEO2{file_ind} = boozer_s_NEO2_tmp(1:data_ctr);
    boozer_psi_pr_NEO2{file_ind} = psi_pr_hat_NEO2_tmp(1:data_ctr).*Bref_NEO2_tmp(1:data_ctr);
    aiota_NEO2{file_ind} = aiota_NEO2_tmp(1:data_ctr);
    bcovar_tht_NEO2{file_ind} = bcovar_tht_NEO2_tmp(1:data_ctr);
    bcovar_phi_NEO2{file_ind} = bcovar_phi_NEO2_tmp(1:data_ctr);
    avb2_NEO2{file_ind} = avbhat2_NEO2_tmp(1:data_ctr).*(Bref_NEO2_tmp(1:data_ctr).^2);
    Er_NEO2{file_ind} = Er_NEO2_tmp(1:data_ctr);
    avEparB_ov_avb2_NEO2{file_ind} = avEparB_ov_avb2_NEO2_tmp(1:data_ctr);
    avnabpsi_NEO2{file_ind} = avnabpsi_NEO2_tmp(1:data_ctr);
    Mt_e_NEO2{file_ind} = MteOvR_NEO2_tmp(1:data_ctr).*R0_NEO2_tmp(1:data_ctr);
    Mt_d_NEO2{file_ind} = MtdOvR_NEO2_tmp(1:data_ctr).*R0_NEO2_tmp(1:data_ctr);
    TphiNA_tot_NEO2{file_ind} = TphiNA_tot_NEO2_tmp(1:data_ctr);
    TphiNA_spec_NEO2{file_ind} = TphiNA_spec_NEO2_tmp(1:data_ctr);
    Gamma_AX_spec_NEO2{file_ind} = Gamma_AX_spec_NEO2_tmp(1:data_ctr);
    Gamma_AX_Ware_spec_NEO2{file_ind} = Gamma_AX_Ware_spec_NEO2_tmp(1:data_ctr);
    Gamma_NA_spec_NEO2{file_ind} = Gamma_NA_spec_NEO2_tmp(1:data_ctr);
    Gamma_NA_Ware_spec_NEO2{file_ind} = Gamma_NA_Ware_spec_NEO2_tmp(1:data_ctr);
    n_spec_NEO2{file_ind} = n_spec_NEO2_tmp(1:data_ctr);
    T_spec_NEO2{file_ind} = T_spec_NEO2_tmp(1:data_ctr);
    Qflux_AX_spec_NEO2{file_ind} = Qflux_AX_spec_NEO2_tmp(1:data_ctr);
    Qflux_AX_Ware_spec_NEO2{file_ind} = Qflux_AX_Ware_spec_NEO2_tmp(1:data_ctr);
    Qflux_NA_spec_NEO2{file_ind} = Qflux_NA_spec_NEO2_tmp(1:data_ctr);
    Qflux_NA_Ware_spec_NEO2{file_ind} = Qflux_NA_Ware_spec_NEO2_tmp(1:data_ctr);
  end

  %% Compute integral NTV torque and surface area

  TphiNA_int_tot_NEO2 = cell(num_data_NEO2,1);
  TphiNA_int_ele_NEO2 = cell(num_data_NEO2,1);
  TphiNA_int_io_NEO2 = cell(num_data_NEO2,1);
  surf_area_NEO2 = cell(num_data_NEO2,1);

  for file_ind = 1:num_data_NEO2
    % input
    boozer_s_NEO2_tmp = boozer_s_NEO2{file_ind};
    boozer_psi_pr_NEO2_tmp = boozer_psi_pr_NEO2{file_ind};
    aiota_NEO2_tmp = aiota_NEO2{file_ind};
    bcovar_tht_NEO2_tmp = bcovar_tht_NEO2{file_ind};
    bcovar_phi_NEO2_tmp = bcovar_phi_NEO2{file_ind};
    avb2_NEO2_tmp = avb2_NEO2{file_ind};
    avnabpsi_NEO2_tmp = avnabpsi_NEO2{file_ind};
    TphiNA_tot_NEO2_tmp = TphiNA_tot_NEO2{file_ind};
    num_radial_pts = numel(TphiNA_tot_NEO2_tmp);
    TphiNA_ele_NEO2_tmp = zeros(num_radial_pts,1);
    TphiNA_io_NEO2_tmp = zeros(num_radial_pts,1);
    for k = 1:num_radial_pts
      TphiNA_ele_NEO2_tmp(k) = TphiNA_spec_NEO2{file_ind}{k}(1);
      TphiNA_io_NEO2_tmp(k) = TphiNA_spec_NEO2{file_ind}{k}(2);
    end
    % local array
    TphiNA_int_tot_NEO2_tmp = zeros(size(TphiNA_tot_NEO2_tmp));
    TphiNA_int_ele_NEO2_tmp = zeros(size(TphiNA_ele_NEO2_tmp));
    TphiNA_int_io_NEO2_tmp = zeros(size(TphiNA_io_NEO2_tmp));
    % compute integral torque
    TphiNA_int_tot_NEO2_tmp(1) = 0;
    TphiNA_int_ele_NEO2_tmp(1) = 0;
    TphiNA_int_io_NEO2_tmp(1) = 0;
    for es_ind = 2:numel(TphiNA_tot_NEO2_tmp)
      TphiNA_int_tot_NEO2_tmp(es_ind) = ...
          trapz(boozer_s_NEO2_tmp(1:es_ind),TphiNA_tot_NEO2_tmp(1:es_ind).*...
          (aiota_NEO2_tmp(1:es_ind).*bcovar_tht_NEO2_tmp(1:es_ind)+bcovar_phi_NEO2_tmp(1:es_ind))./...
          avb2_NEO2_tmp(1:es_ind));
      TphiNA_int_ele_NEO2_tmp(es_ind) = ...
          trapz(boozer_s_NEO2_tmp(1:es_ind),TphiNA_ele_NEO2_tmp(1:es_ind).*...
          (aiota_NEO2_tmp(1:es_ind).*bcovar_tht_NEO2_tmp(1:es_ind)+bcovar_phi_NEO2_tmp(1:es_ind))./...
          avb2_NEO2_tmp(1:es_ind));
      TphiNA_int_io_NEO2_tmp(es_ind) = ...
          trapz(boozer_s_NEO2_tmp(1:es_ind),TphiNA_io_NEO2_tmp(1:es_ind).*...
          (aiota_NEO2_tmp(1:es_ind).*bcovar_tht_NEO2_tmp(1:es_ind)+bcovar_phi_NEO2_tmp(1:es_ind))./...
          avb2_NEO2_tmp(1:es_ind));
    end
    TphiNA_int_tot_NEO2_tmp = (4*pi^2)*boozer_psi_pr_NEO2_tmp.*TphiNA_int_tot_NEO2_tmp;
    TphiNA_int_ele_NEO2_tmp = (4*pi^2)*boozer_psi_pr_NEO2_tmp.*TphiNA_int_ele_NEO2_tmp;
    TphiNA_int_io_NEO2_tmp = (4*pi^2)*boozer_psi_pr_NEO2_tmp.*TphiNA_int_io_NEO2_tmp;
    % store data for plotting
    TphiNA_int_tot_NEO2{file_ind} = TphiNA_int_tot_NEO2_tmp;
    TphiNA_int_ele_NEO2{file_ind} = TphiNA_int_ele_NEO2_tmp;
    TphiNA_int_io_NEO2{file_ind} = TphiNA_int_io_NEO2_tmp;
    surf_area_NEO2{file_ind} = (4*(pi^2)*avnabpsi_NEO2_tmp./avb2_NEO2_tmp).*...
        abs(boozer_psi_pr_NEO2_tmp.*(aiota_NEO2_tmp.*bcovar_tht_NEO2_tmp+bcovar_phi_NEO2_tmp));
  end

end
%% Data for Sergei

function plot_2spec_export(boozer_s, TphiNA_tot, TphiNA_int_tot, TphiNA_int_ele, TphiNA_int_io)
  ind_data=1;
  disp(['integral NTV torque = ',num2str(1e-7*TphiNA_int_tot{ind_data}(end)),' Nm'])
  disp(['integral NTV torque (ele) = ',num2str(1e-7*TphiNA_int_ele{ind_data}(end)),' Nm'])
  disp(['integral NTV torque (io) = ',num2str(1e-7*TphiNA_int_io{ind_data}(end)),' Nm'])


  fid_vphiref = fopen('vphiref.in');
  vphiref_str = fgetl(fid_vphiref);
  vphiref_val_str = regexprep(vphiref_str,'[^0-9.\-]','');
  vphiref = str2num(vphiref_val_str);

  NTV_torque_int_output=[vphiref,1e-7*TphiNA_int_tot{ind_data}(end),1e-7*TphiNA_int_ele{ind_data}(end),1e-7*TphiNA_int_io{ind_data}(end)];
  save('NTV_torque_int.dat','NTV_torque_int_output','-ascii')
end
