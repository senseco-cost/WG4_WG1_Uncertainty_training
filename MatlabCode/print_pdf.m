function print_pdf(figure_handle, path, filename)

        if nargin == 2
            fullfilename = path;
        else
            fullfilename = [path filename '.pdf'];
        end

        %fillPage(figure_handle, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
        fillPage(figure_handle, 'margins', [0 0 0 0], 'papersize', [30 20]);
        print(figure_handle, '-dpdf', '-r600', fullfilename)
end