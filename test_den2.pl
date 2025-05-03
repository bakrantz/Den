use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Den";
use Den;
my $d = Den -> new() 
            -> open_parameter_file('den_parameter_file_2.txt')  
            -> calc_den() 
            -> extract_wavelength_signal()
            -> save_data('spectrum_MELT.txt');
