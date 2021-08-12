function convert_to_TRAVISformat(file_in, file_out)
  %~ file_in  = '/proj/plasma/Neo2/Interface/Profiles/w7x-m111-b3-i1/prf.txt';
  %~ file_out = '/proj/plasma/Neo2/Interface/Profiles/w7x-m111-b3-i1/profiles_for_travis.dat';

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

  skip = 1;
  %% open output file and write first lines
  fid = fopen(file_out,'w');
  fprintf(fid, 'CC                r/a               ne,m^3               Te,keV             Zeff \n');
  for j = 1:skip
    %fprintf(fid, '%s\n', data_c{j});
    fprintf(fid, 'numberOfPoints %d\n', k - skip);
  end

  data_c = data_c(1+skip:end);
  data_m = cellfun(@(x) str2num(x), data_c, 'UniformOutput', false).';
  data_m = cell2mat(data_m);

  cols = [1,2,4,7];   % r/a, ne, Te, Zeff

  data_m = data_m(:, cols);

  % Conversion
  data_m(:,2) = data_m(:,2) * 1e+20;
  %data_m(:,3) = data_m(:,3);

  %cols_sort = [1, 2, 4, 3];
  %data_m = data_m(:, cols_sort);

  for j = 1:size(data_m,1)
    fprintf(fid,'%20.10E %20.10E %20.10E %20.10E\n', data_m(j,1), data_m(j,2), data_m(j,3), data_m(j,4));
  end
  fclose(fid);
end
