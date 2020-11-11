function print_pdf(figure_handle, path, filename)

    figure_handle.Color = 'white';
    export_fig(figure_handle, [path filename '.pdf']);

end

%        fillPage(figure_handle, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
%        print(figure_handle, '-dpdf', '-r600', [path filename '.pdf'])