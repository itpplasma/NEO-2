%
% input:
% ------
% file_base: base name of input file (without file extension), and for
%   output files. Output file names will be 'file_base-n#mode[.file_ext]',
%   where '#mode' is replaced with the mode number.
% file_ext: file extension of input file. If not passed or empty, it is
%   assumed there is no file extension.
% file_descr: string, additional information about fily type to read.
%   Defaults to 'cos_harm'.
function extract_pert_field(file_base, file_ext, file_descr)
  if nargin < 3 || isempty(file_descr)
    file_descr = 'cos_harm'; % additional info (some extra info given by the filename)
  end

  %file_base = 'aug_2_n0_shortto_m12-pert-2-3'%'aug_2_rmp-n0.bc'; % 'tok-synch2.bc'; % data file base name (identifier)
  %file_ext = 'bc'; %extension

  % I try to read the file
  if nargin < 2 || isempty(file_ext)
    file_in = file_base;
    file_ext = []; % in case it is not set.
  else
    file_in = [file_base, '.', file_ext];
  end
  fid = fopen(file_in);

  k = 1;
  tline = fgetl(fid);
  data_c{k} = tline;
  while ischar(tline)
    % disp(tline)
    tline = fgetl(fid);
    if ischar(tline)
      k = k + 1;
      data_c{k} = tline;
    end
  end
  fclose(fid);
  % all lines are now in data_c

  % find m0b and n0b
  for k = 1:numel(data_c)
    tline = data_c{k};
    if strfind(tline,'m0b')
      k_def = k + 1;
      break
    end
  end

  % extract m0b and n0b
  tline = data_c{k_def};
  dline = str2num(tline);
  m0b = dline(1);
  n0b = dline(2);
  nsurf = dline(3);
  nper = dline(4);

  % check the first block of the spectrum (negative n- or m-values)
  % (NEO-2 needs B field spectrum of the from n>=0 && m>-inf)
  spec = data_c(k_def+5:k_def+5+(n0b+1+m0b*(2*n0b+1))-1);
  line_conv = str2num(spec{1});
  spec_num = zeros(n0b+1+m0b*(2*n0b+1),size(line_conv,2));
  spec_num(1,:) = line_conv;
  for k=2:(n0b+1+m0b*(2*n0b+1))
    line_conv = str2num(spec{k});
    spec_num(k,:) = line_conv;
  end
  %get the highest toroidal perturbation mode number
  n_pert_max = max(unique(spec_num(:,2)));
  %check m>-inf
  check_m_symm = sum(flipud(unique(spec_num(:,1)))+unique(spec_num(:,1)))==0;
  %check n>=0
  check_n_symm = sum(flipud(unique(spec_num(:,2)))+unique(spec_num(:,2))-n0b)==0;
  %specify switch for mapping Foruier coefficients
  if check_m_symm && check_n_symm % this is the wanted spectrum
    sw_spec = 1;
  else
    %check m>=0
    check_m_symm = sum(flipud(unique(spec_num(:,1)))+unique(spec_num(:,1))-m0b)==0;
    %check n>-inf
    check_n_symm = sum(flipud(unique(spec_num(:,2)))+unique(spec_num(:,2)))==0;
    if check_m_symm && check_n_symm % Fourier coefficients must be mapped
      sw_spec = 2;
    else
      error('Check your Boozer file (mode numbers)!')
    end
  end
  % change now n0b and m0b according to perturbation
  dline(2) = 0; % only 1 toroidal mode number per file
  switch sw_spec
  case 1
    dline(1) = m0b;
  case 2
    % Requires additional header for the n equal zero mode.
    data_c0(1:k_def)=data_c(1:k_def);
    data_c0{k_def} = sprintf(' %d    %d    %d   %d   %.6e   %.5f   %.5f', dline);
    dline(1) = m0b * 2; % because of negative perturbations
  end
  data_c{k_def} = sprintf(' %d    %d    %d   %d   %.6e   %.5f   %.5f', dline);

  % open the output files for the different
  % (toroidal) perturbation mode numbers
  file_names=cell(1,n_pert_max+1);
  fids=cell(1,n_pert_max+1);
  for n_pert=0:n_pert_max
    file_out = [file_base,'-n',num2str(n_pert)];
    if ~isempty(file_ext)
      file_out=[file_out,'.',file_ext]; %#ok<AGROW>
    end
    file_names{n_pert+1}=file_out;
    fids{n_pert+1}=fopen(file_out,'w');
    if n_pert==0
      for kl = 1:k_def-1
        fprintf(fids{n_pert+1}, '%s\n', data_c0{kl});
      end
      line_conv=str2num(data_c0{k_def});
      line_conv(4) = 1;
      fprintf(fids{n_pert+1},'  %2d    %2d   %4d    %2d    %9.8f     %9.8f     %9.8f\n', line_conv);
    else
      for kl = 1:k_def
        fprintf(fids{n_pert+1}, '%s\n', data_c{kl});
      end
    end
  end

  % separate spectra accodring to toroidal mode number and
  % write data into separate files
  k_start = k_def;
  for n = 1:nsurf
    head = data_c(k_start+1:k_start+4);

    %keyboard
    if strcmp(file_descr, 'cos_harm')
      head{end} =  ['    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]'];
    end
    spec = data_c(k_start+5:k_start+5 + (n0b+1+m0b*(2*n0b+1)) - 1);


    % write the original header for flux surface
    for n_pert=0:n_pert_max
      if n_pert==0
        for kl = 1:numel(head)-2
          fprintf(fids{n_pert+1},'%s\n',head{kl});
        end
        line_conv = str2num(head{numel(head)-1});
        line_conv(3) = line_conv(3)*nper;
        line_conv(6) = line_conv(6)*nper;
        fprintf(fids{n_pert+1}, '   %9.8e   %9.8e   %9.8e   %9.8e  %9.8e  %9.8e\n', line_conv);
        fprintf(fids{n_pert+1}, '%s\n', head{numel(head)});
      else
        for kl = 1:numel(head)
          fprintf(fids{n_pert+1}, '%s\n', head{kl});
        end
      end
    end

    %convert spectrum to a matrix
    line_conv=str2num(spec{1});
    ncol_spec=size(line_conv,2);
    spec_num=zeros((n0b+1+m0b*(2*n0b+1)),ncol_spec);
    spec_num(1,:)=line_conv;
    for k=2:(n0b+1+m0b*(2*n0b+1))
      line_conv=str2num(spec{k});
      spec_num(k,:)=line_conv;
    end

    %extract single (toroidal) perturbations
    switch sw_spec
    case 1 % this is the wanted form of the spectrum
    case 2 % here Fourier coefficients must be mapped
      for n_pert=0:n_pert_max
        if n_pert == 0 % this is the axisymmetric field
          L = spec_num(:,2) == 0;
          data_6cols = spec_num(L,:);
          switch size(data_6cols,2)
          case 6
            if strcmp(file_descr,'cos_harm')
              empty_col = zeros(size(data_6cols(:,1)));
              data = [data_6cols(:,1), data_6cols(:,2),  data_6cols(:,3), empty_col, empty_col, data_6cols(:,4),  empty_col, data_6cols(:,5), data_6cols(:,6), empty_col];
            end
          case 10
            data=spec_num(L,:);
          otherwise
            error('Invalid format in input .bc file.')
          end


          data(:,7)=data(:,7)/nper;
          data(:,8)=data(:,8)/nper;

          for kl = 1:size(data,1)
            fprintf(fids{n_pert+1},'   %2d   %2d  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e\n', data(kl,:));
          end
        else % perturbation field - coefficients must be mapped
          try
            data_new=zeros(2*m0b+1,ncol_spec);
            L1=spec_num(:,2)==-n_pert;
            L_m0=spec_num(:,1)==0;

            data_new(1:m0b,:)=[-spec_num(~L_m0&L1,1),...
                abs(spec_num(~L_m0&L1,2)),spec_num(~L_m0&L1,3:end)];
            %sin()-contribution
            data_new(1:m0b,4:2:end)=-data_new(1:m0b,4:2:end);
            L2=spec_num(:,2)==n_pert;

            data_new(m0b+1:2*m0b+1,:)=[spec_num(L2,1),...
                abs(spec_num(L2,2)),spec_num(L2,3:end)];
            [~,ind]=sort(data_new(:,1));
            data_new=data_new(ind,:);

            switch size(data_new,2)
            case 6
              if strcmp(file_descr,'cos_harm')
                empty_col = zeros(size(data_new(:,1)));
                data_new_col10 = [data_new(:,1), data_new(:,2), data_new(:,3), empty_col, empty_col, data_new(:,4), empty_col, data_new(:,5), data_new(:,6), empty_col];
              end
            case 10
              data_new_col10 = data_new;
            otherwise
              error('Invalid format in input .bc file.')
            end

            for kl = 1:size(data_new,1)
              fprintf(fids{n_pert+1},'  %3d   %2d  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e  %15.8e\n', data_new_col10(kl,:));
            end
          end
        end
      end
    end

    if strcmp(file_descr,'cos_harm')
      k_start = k_start+5 + ((m0b+1)*(2*n0b+1)) - 1;
    else
      k_start = k_start+5 + (n0b+1+m0b*(2*n0b+1)) - 1;
    end

  end

  % close the output files for the different
  % (toroidal) perturbation mode numbers
  for n_pert=0:n_pert_max
    fclose(fids{n_pert+1});
  end
end
