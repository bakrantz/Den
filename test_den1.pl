use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Den"; #Make a path to the module
use Den;
my $d = Den -> new() 
            -> open_parameter_file('den_parameter_file_1.txt') 
            -> open_raw_data('01.PRN') 
            -> calc_den() 
            -> extract_signal()
            -> save_data('01_MELT.txt');

