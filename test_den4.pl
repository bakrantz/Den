use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Den"; #Make a path to the module
use Den;
my $den = Den -> new() 
              -> open_parameter_file('den_parameter_file_4.txt', chr(44)) #set to open comma delimited parm file
              -> calc_den()
              -> extract_wavelength_signal('FRET', 'txt', chr(9)) #set to extract from tab delimited spectra txt files
              -> save_data('spectrum_FRET_2_MELT.txt');

