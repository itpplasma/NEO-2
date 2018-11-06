function s = h52struct(fname,finfo,fbase,s,h5id)

if nargin < 2 || isempty(finfo), finfo = hdf5info(fname); end
if nargin < 3 || isempty(fbase), fbase = 'finfo.GroupHierarchy'; end
if nargin < 4 || isempty(s),     s     = struct(); end
if nargin < 5 || isempty(h5id),  h5id  = H5F.open(fname); end

fgroups = [fbase,'.Groups'];
fgroupsname = ['{',fgroups,'.Name','}'];

if ~isempty(eval(fgroups))
    groups = eval(fgroupsname);
    for kg = 1:numel(groups)
        field = groups{kg}(2:end);
        field = regexprep(field,'[^\w]','_');
        fbase_n=[fgroups,sprintf('(%i)',kg)];
        s.(field) = struct();
        s.(field) = h52struct(fname,finfo,fbase_n,s.(field),h5id);
    end
end

fdatasets = [fbase,'.Datasets'];
fdataname = ['{',fdatasets,'.Name','}'];
if ~isempty(eval(fdatasets))
    datasets = eval(fdataname);
    for kd = 1:numel(datasets)
        dataset = datasets{kd};
        d = strsplit(dataset(2:end),filesep);
        d = regexprep(d,'[^\w]','_');
        
        % High-Level
        %data = h5read(fname,dataset);
        
        % Low-Level
        dset_id = H5D.open(h5id,dataset);
        data = H5D.read(dset_id);
        H5D.close(dset_id);        
        
        if isa(data,'hdf5.h5string')
            s.(d{end}) = data.Data;
        else
            s.(d{end}) = data;
        end
    end
end