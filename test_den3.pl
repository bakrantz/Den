use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Den"; #Make a path to the module
use Den;
my $den = Den -> new() 
              -> open_parameter_file('den_parameter_file_3.txt') #open tab delimited parm file by default
              -> calc_den()
              -> extract_wavelength_signal()
              -> save_data('spectrum_FRET_MELT.txt');


