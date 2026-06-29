function [data] = read_file_into_cell_array(filename)
  fid = fopen(filename);

  k = 1;
  tline = fgetl(fid);
  data{k} = tline;
  while ischar(tline)
    tline = fgetl(fid);
    if ischar(tline)
      k = k + 1;
      data{k} = tline;
    end
  end

  fclose(fid);
end
