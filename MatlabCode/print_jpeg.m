function print_jpeg(figure_handle, path, filename)
%         fillPage(figure_handle,  'margins', [0 0 0 0], 'papersize', [11 8.5]);
%         print(figure_handle,  '-djpeg', '-r600', [path filename '.jpg'])
        
    figure_handle.Color = 'white';
    export_fig(figure_handle, [path filename '.jpg'], '-jpg', '-r600');
        
end