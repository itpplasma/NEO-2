function export_plot( name, varargin )

    plot_device = '-depsc2';
    ext = '.eps';
    plot_name = [ name, ext];
    print(plot_device, plot_name);
    convert_eps_2_pdf(plot_name);
    
end