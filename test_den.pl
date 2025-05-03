use lib "/Users/bakrantz/Documents/perl/lib/perl/5.8/Den"; #Make a path to the module
use Den;
my $d = Den -> new() 
            -> open_parameter_file('den_parameter_file.txt') 
            -> conjoin_files('072808A.PRN','072808B.PRN','conjoined.txt') 
            -> calc_den() 
            -> extract_signal()
            -> save_data('072808_MELT.txt');
