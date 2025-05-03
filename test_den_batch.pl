use Den;
foreach my $prefix ('03','02','01') {
my $den_obj = Den 	-> new()
					-> open_parameter_file('den_parameter_file.txt')
					-> open_raw_data($prefix.'.PRN')
					-> calc_den()
					-> extract_signal()
					-> save_data($prefix.'_MELT.txt')
};