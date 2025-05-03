#####################################################################################
#################                   Den.pm                           ################
#################             Krantz Lab       July 2008             ################
#################                              Rev. Dec 2022         ################
#####################################################################################
################# Module to process denaturant or ligand titration data #############
#####################################################################################
package Den;
################################# THE CONSTRUCTOR ###################################
sub new {
  my $class                                   = shift;
  my $self                                    = {};                 #the anon hash
     $self->{DATA}                            = [];
     $self->{RAW_DATA}                        = []; 
     $self->{RAW_SPECTRUM}                    = [];
     $self->{DEN}                             = [];
     $self->{LIGAND}                          = [];
     $self->{SIGNAL}                          = [];
     $self->{TITRATION}                       = [];
     $self->{LIGAND_TITRATION}                = [];
     $self->{LIGAND_STOCK_CONCENTRATIONS}     = [];
     $self->{DEN_NAME}                        = undef;    
     $self->{LIGAND_NAME}                     = undef;
     $self->{SIGNAL_NAME}                     = undef;
     $self->{HEADER}                          = undef;
     $self->{CUVETTE_VOLUME}                  = undef;
     $self->{TITRANT_CONC}                    = undef;
     $self->{CUVETTE_LIGAND_CONCENTRATION}    = undef;
     $self->{CUVETTE_CONC}                    = undef;
     $self->{TITRATION_PROGRAM}               = undef;
     $self->{FLUORESCENCE_MODE}               = undef;
     $self->{TITRATION_MODE}                  = undef;
     $self->{TIME_OFFSET}                     = undef;
     $self->{TIME_INIT}                       = undef;
     $self->{TIME_PER_BIN}                    = undef;
     $self->{TIME_START}                      = undef;
     $self->{TIME_END}                        = undef;
     $self->{WAVELENGTH}                      = undef; 
     $self->{WAVELENGTH_DONOR}                = undef; 
     $self->{WAVELENGTH_ACCEPTOR}             = undef; 
     $self->{WAVELENGTH_WIDTH}                = undef;
     $self->{PARAMETER_FILE}                  = undef;
     $self->{FILE_OUT}                        = undef;
     $self->{FILE}                            = undef;
     $self->{FILE_SPECTRUM}                   = undef;
     $self->{FILE_SPECTRUM_EXTENSION}         = undef;
     $self->{FILE_SPECTRUM_HEADER_LENGTH}     = undef;
     $self->{DELIMITER}                       = undef; 
     $self->{VERBOSE}                         = undef;
  return bless($self, $class)                                  #return thy self
}

############################### METHOD SUBS ########################################

#Scalar properties can be loaded into the module from a parameter file
# Generally, the file is in KEY tab value format for the parameters but you may set another delimiter
# (other than tab, which is default)
sub open_parameter_file {
  my ($self, $file, $delimiter) = @_;
  $file = parameter_file($self) unless $file;
  $delimiter = delimiter($self) unless $delimiter;
  $delimiter = delimiter($self, chr(9)) unless $delimiter;
  delimiter($self, $delimiter);
  open(FH, '<', parameter_file($self, $file)) || die "Cannot open file called $file in open_parameter_file method. $!";
  while (my $line = <FH>) {
    my ($key, $value) = (undef, undef);
    ($key, $value) = split(/\s*$delimiter\s*/, ts($line)) if ($delimiter ne chr(9)); #split on delimiter between keys and values
    ($key, $value) = split(/\s+/, ts($line)) if ($delimiter eq chr(9)); #sloppy split on tab with whitespace regexp
    $self -> {$key} = $value;
  };
  close(FH);
  return $self
}

#TITRATION_PROGRAM sets the titration protocol used, for example, two segments are
#as follows in its shorthand form, i.e., 30_25_48_64.
#The first segment has 30 pts of 25 uL titrant volume; the second segment has 48 pts of 64 uL titrant volume.
#In a table it would be:
#
#SEG. 	PTS.	VOLUME (UL)
#============================
#SEG1	30	25
#SEG2	48	64
#============================
#
#NOTE the TITRATION_PROGRAM must separate the values for the segments using only an underscore
#and no spaces or other characters—only numbers and underscores.
#A simple one-segment program would be 153_25, which is 153 points of 25 uL titrant volume per point
sub init_titration_program {
  my ($self, $titration_program) = @_;
  my $us = '_'; #numerical values in titration program string are always separated by underscores 
  $titration_program = titration_program($self) unless $titration_program;
  if ($titration_program) {
    titration_program($self, $titration_program);
    print "titration_program is $titration_program\n" if verbose($self);
  }
  else { die "Aborting init_titration_program because no titration_program is available. Please set in parameter file.\n" };
  my @titr_prog = split(/$us/, $titration_program);
  my $segments = int(scalar(@titr_prog)/2) - 1;
  if ($segments < 0) { die "Aborting init_titration_program because $titration_program is not in correct steps-underscore-volume format."}
  my @pts = (); #number of points in each segment 
  my @vols = (); #volume titrated per point in each segment
  my @titration = (); #expanded single-wide array of length total titration points containing the value of the respective titration volumes 
  for (0..$segments) {
    my $pp = shift(@titr_prog);
    print "Segment $_ has $pp titrations\n"  if verbose($self);
    push @pts, $pp;
    my $vv = shift(@titr_prog);
    print "Segment $_ dispenses a volume of $vv uL per step\n"  if verbose($self);
    push @vols, $vv; 
  };
  for my $i (0..$segments) { if ($pts[$i]) { for (1..$pts[$i]) { push @titration, $vols[$i] } } };
  titration($self, @titration);
  return $self
}

#Calculates the denaturant concentration array for a constant cuvette volume titration
#A set volume is removed from the cuvette and the same volume of titrant with denaturant is added to the cuvette
# where protein concentration remains constant and cuvette is under continuous mixing
sub calc_den {
  my ($self, $titrant_conc, $cuvette_conc, $cuvette_volume, $titration_program) = @_;
  my $titration = []; #ref to array for the array of titration volumes given in the titration_program
  $titrant_conc      = titrant_conc($self)      unless $titrant_conc; #titrant [Denaturant] in Molar
  $cuvette_conc      = cuvette_conc($self)      unless $cuvette_conc; #cuvette [Denaturant] in Molar
  $cuvette_volume    = cuvette_volume($self)    unless $cuvette_volume; #cuvette volume in uL
  $titration_program = titration_program($self) unless $titration_program; 
  if ($titration_program) {
    init_titration_program($self, $titration_program);
    $titration      = [ titration($self) ];
    titration_program($self, $titration_program);
  }
  else { die "Aborting no titration_program available. Please set in parameter file.\n" };
  $titrant_conc     = titrant_conc($self, $titrant_conc);
  print "titrant_conc is $titrant_conc M\n"  if verbose($self);
  $cuvette_conc = cuvette_conc($self, $cuvette_conc);
  print "cuvette_conc is $cuvette_conc M\n"  if verbose($self);
  $cuvette_volume   = cuvette_volume($self, $cuvette_volume);
  print "cuvette_volume is $cuvette_volume uL\n"  if verbose($self);
  $titration        = [titration($self, @$titration)];
  my @steps = @$titration;
  my $number_of_steps = @steps;
  my @den = ();
  push @den, $cuvette_conc; #this is the initial point (cuvette conc.) before any titration has started
  print "Number of titration steps: $number_of_steps\n\n"  if verbose($self);
  my $dn = den_name($self) if den_name($self);
  $dn = den_name($self, 'DEN') unless $dn;
  print "Point\t".'['.$dn.']'."\n"  if verbose($self);
  print "====================\n"  if verbose($self);
  print "0\t".$cuvette_conc."\n"  if verbose($self);  
  for my $i (0..$#steps) {
    # Calculate the titration
    $cuvette_conc = ($cuvette_conc * ($cuvette_volume - $steps[$i]) + $titrant_conc * $steps[$i]) / $cuvette_volume;
    push @den, $cuvette_conc;
    my $titr = $i + 1;
    print "$titr\t$cuvette_conc\n"  if verbose($self);  
  };
  print "\n" if verbose($self);
  titration_mode($self, 'DENATURANT'); #mode is set to 'DENATURANT' so that stitch_data and save_data know what the x-column is
  den($self, @den);
  return $self
}

#For a non-denaturant ligand titration, use calc_ligand to generate the ligand concentrations for each titration point
# $ligand_stock_concentrations is a reference to an array of ligand stock concentrations in the order used in the titration
# $ligand_titration is a reference to an array of arrays of the volumes titrated for each stock concentration
# a $ligand_titration example is [ [1, 2.5, 6], [1, 2.5, 6], [1, 2.5, 6], [1, 2.5, 6] ]
# which delivers 1, 2.5 and 6 microliter points per stock concentration (of which there are four)
# $cuvette_volume is the volume of the cuvette in microliters.
# $cuvette_ligand_concentration begins as the starting concentration of ligand
# in the cuvette. This cuvette concentration increases as ligand is titrated during the experiment
# $ligand_name is the name of the ligand used; it labels the appropriate data column in the output file in save_data
sub calc_ligand {
  my ($self, $ligand_stock_concentrations, $ligand_titration, $cuvette_volume, $cuvette_ligand_concentration, $ligand_name) = @_;
  my @ligand = ();
  my $point = 0;
  $ligand_stock_concentrations = [ ligand_stock_concentrations($self) ] unless $ligand_stock_concentrations;
  $ligand_titration = [ ligand_titration($self) ] unless $ligand_titration;
  $cuvette_volume = cuvette_volume($self) unless $cuvette_volume;
  $cuvette_ligand_concentration = cuvette_ligand_concentration($self) unless $cuvette_ligand_concentration;
  $ligand_name = ligand_name($self) unless $ligand_name;
  $ligand_name = 'LIGAND' unless $ligand_name;  
  push @ligand, $cuvette_ligand_concentration; #this is the initial point prior to any additions
  print "Point\t".'['.$ligand_name.']'."\n"           if verbose($self);
  print "====================\n"                      if verbose($self);
  print "$point\t".$cuvette_ligand_concentration."\n" if verbose($self);  
  for (my $i = 0; $i < scalar(@$ligand_stock_concentrations); $i++) {
    my $stock_conc = $ligand_stock_concentrations -> [$i];
    my $volumes = $ligand_titration -> [$i];
    for (my $j = 0; $j < scalar(@$volumes); $j++) {
      my $volume = $volumes -> [$j];
      my $moles_ligand_cuvette = $cuvette_volume * $cuvette_ligand_concentration;
      my $moles_ligand_added = $volume * $stock_conc;
      $cuvette_volume += $volume;
      $cuvette_ligand_concentration = ($moles_ligand_cuvette + $moles_ligand_added) / $cuvette_volume;
      push @ligand, $cuvette_ligand_concentration;
      $point++;
      print "$point\t".$cuvette_ligand_concentration."\n" if verbose($self);
    };   
  };
  ligand_stock_concentrations($self, @$ligand_stock_concentrations);
  ligand_titration($self, @$ligand_titration);
  cuvette_volume($self, $cuvette_volume);
  cuvette_ligand_concentration($self, $cuvette_ligand_concentration);
  ligand_name($self, $ligand_name);
  #important to set titration_mode in module so stitch_data and save_data methods
  # know its a ligand titration and not a denaturant titration
  titration_mode($self, 'LIGAND'); 
  ligand($self, @ligand);
  return $self
}

#For a timed denaturant titration, the time versus signal file is opened as raw data
# and stored in the module under raw_data
sub open_raw_data { 
  my ($self, $file, $del) = @_;
  $file = file($self) unless $file;
  $del = delimiter($self) unless $del;
  $del = delimiter($self, chr(9)) if !$del;
  delimiter($self, $del);
  my @raw_data = ();
  open(FH, '<', file($self, $file)) || die "Cannot open $file in open_raw_data. $!";
  print "Opening raw timecourse data file $file\n" if verbose($self);
  while (my $line = <FH>) {
    my ($time, $sig) = split(/\s*$del\s*/, ts($line));
    push @raw_data, [$time, $sig];
  };
  close(FH);
  raw_data($self, @raw_data);
  return $self
}

#When the individual denaturant titration points are spectra, this method opens each spectrum
# and saves the data in raw_spectrum. The $header_length parameter is the number of lines at the
# beginning of the spectrum file that precedes the actual wavlength and signal data
sub open_spectrum {
  my ($self, $file_spectrum, $del, $header_length) = @_;
  $file_spectrum = file_spectrum($self) unless $file_spectrum;
  $del = delimiter($self) unless $del;
  $del = delimiter($self, chr(9)) if !$del;
  delimiter($self, $del);
  $header_length = file_spectrum_header_length($self) unless $header_length; 
  my @raw_spectrum = ();
  open(FH, '<', file_spectrum($self, $file_spectrum)) || die "Cannot open $file_spectrum in open_spectrum. $!";
  print "Opening spectrum file $file_spectrum\n" if verbose($self);
  print "Header length of spectrum file is $header_length lines\n" if verbose($self);
  if ($header_length) {
    foreach (1..$header_length) {
      my $line = <FH>;
      $line = ts($line);
      print "Line $_ of header is $line\n" if verbose($self);
    };
  };
  while (my $line = <FH>) {
    my ($wavelength, $signal) = split(/\s*$del\s*/, ts($line));
    push @raw_spectrum, [$wavelength, $signal];
    print "$wavelength\t$signal\n" if verbose($self);
  };
  close(FH);
  raw_spectrum($self, @raw_spectrum);
  return $self
}

#Joins together two timed denaturant titration raw data files,
# where the second file's times are offset by time_offset property
sub conjoin_files {
  my ($self, $file_a, $file_b, $conjoined_file, $time_offset, $delimiter) = @_;
  my @file = ($file_a, $file_b);
  $time_offset = time_offset($self) unless $time_offset;
  time_offset($self, $time_offset);
  $conjoined_file = file_out($self) unless $conjoined_file;
  $conjoined_file = file_out($self, "conjoined_$file_a".'_'."$file_b".'.txt') unless $conjoined_file;
  $delimiter      = delimiter($self) unless $delimiter;
  $delimiter      = delimiter($self, chr(9)) if !$delimiter;
  delimiter($self, $delimiter);
  my @raw_data = ();
  my @aoas = ();
  for (0..1) { push @aoas, [ file_to_aoa(file($self, $file[$_]), $delimiter) ] }; 
  #offset aoa[1] and then join together
  foreach my $i (@{ $aoas[1] }) { $i->[0] += $time_offset };
  foreach my $j (0..1) { foreach my $i (@{ $aoas[$j] }) { push @raw_data, $i } };
  @raw_data = sort { $a->[0] <=> $b->[0] } @raw_data;
  raw_data($self, @raw_data);
  save_data($self, $conjoined_file, $delimiter, \@raw_data);
  return $self
}

#Opens a file and returns an array of arrays (aoa)
sub file_to_aoa {
  my ($file, $delimiter) = @_;
  $file = file($self) unless $file;
  $delimiter = delimiter($self) unless $delimiter;
  $delimiter = delimiter($self, chr(9)) if !$delimiter;
  delimiter($self, $delimiter);
  my @aoa = ();
  open(FH, '<', file($self, $file)) || die "Cannot open $file in file_to_aoa. $!";
  print "Opening $file to convert to array of arrays prior to conjoining.\n" if verbose($self);
  while (my $line = <FH>) {
    $line = ts($line);
    if ($line =~ /$delimiter/) {
      my @a = split(/\s*$delimiter\s*/, $line);
      push @aoa, [ @a ];
    };
  };
  close(FH);
  return @aoa
}

#Extracts the signal per titration point from a timed denaturant titration
#takes no arguments; all parameters used in the method are set in the parameter file 
sub extract_signal {
  my $self = shift;
  #sort the raw_data
  my @raw_data = sort { $a->[0] <=> $b->[0] } @{ [raw_data($self)] };
  my @signal = ();
  print "Retrieved raw_data array with $#raw_data entries\n"  if verbose($self);
  #get arguments
  my $time_init = time_init($self);
  my $time_per_bin = time_per_bin($self);
  my $time_start = time_start($self);
  my $time_end = time_end($self);
  $time_init = time_init($self, $time_init);
  print "time_init is $time_init\n"        if verbose($self);
  $time_per_bin = time_per_bin($self, $time_per_bin);
  print "time_per_bin is $time_per_bin\n"  if verbose($self);
  $time_start = time_start($self, $time_start);
  $time_end = time_end($self, $time_end);
  #find start and end times if not given
  $time_start = time_start($self, $raw_data[0][0]) unless $time_start;
  $time_end   = time_end($self, $raw_data[-1][0])  unless $time_end;
  print "time_start is $time_start\n"       if verbose($self);
  print "time_end is $time_end\n"           if verbose($self);
  #parse and average raw signal to time bins <= and =>
  my ($this_time, $next_time, $avgsig, $cnt) = ($time_start, $time_start + $time_init, 0, 0);
  for my $i (0..$#raw_data) {
    my ($time, $sig) = ($raw_data[$i][0], $raw_data[$i][1]);
    my $bool = ((($next_time <=> $time ) + 1) * (($time <=> $this_time) + 1));
    if ($bool > 0) {
      $avgsig += $sig;
      $cnt++;
    } 
    else {
      $this_time += $time_per_bin;
      $next_time += $time_per_bin;
      $avgsig = $avgsig / $cnt; #/#division
      push @signal, $avgsig;
      print "Count per bin is $cnt\tAvg signal is $avgsig\n"  if verbose($self);
      ($cnt, $avgsig) = (1, $sig);  #reset cnt and avgsig to next bin
    };
  };
  signal($self, @signal);
  stitch_data($self);
  return $self
}

#Extracts single wavelength fluorescence signal from individual spectrum files for each titration point
#where the spectrum files are indexed 0.txt, 1.txt, 2.txt, 3.txt ... with a 'txt' file extension
#There are two fluorescence modes 'SINGLE_WAVELENGTH' and 'FRET' the wavelength parameters are set in the parameter file
sub extract_wavelength_signal {
  my ($self, $fluorescence_mode, $ext, $delimiter) = @_;
  if ($fluorescence_mode) { fluorescence_mode($self, $fluorescence_mode) } else { $fluorescence_mode = fluorescence_mode($self) }; 
  $ext = file_spectrum_extension($self) unless $ext; #Extension for all the spectrum files
  print "file extension for spectra is $ext\n"             if verbose($self);
  $delimiter = delimiter($self) unless $delimiter;
  $delimiter = chr(9) unless $delimiter;
  print "fluorescence_mode is $fluorescence_mode\n"        if verbose($self);
  my $wavelength = wavelength($self);
  print "wavelength is $wavelength nm\n"                   if (verbose($self) && $wavelength);
  my ($wavelength_width, $wavelength_donor, $wavelength_acceptor) = (wavelength_width($self), wavelength_donor($self), wavelength_acceptor($self));
  print "wavelength_width is $wavelength_width nm\n"       if verbose($self);
  print "wavelength_donor is $wavelength_donor nm\n"       if (verbose($self) && $wavelength_donor);
  print "wavelength_acceptor is $wavelength_acceptor nm\n" if (verbose($self) && $wavelength_acceptor);
  my @signal = ();
  #start and end for averaging window are plus or minus half the width
  my ($start, $end) = ($wavelength - $wavelength_width/2, $wavelength + $wavelength_width/2);
  my ($donor_start, $donor_end) =  ($wavelength_donor - $wavelength_width/2, $wavelength_donor + $wavelength_width/2);
  my ($acceptor_start, $acceptor_end) =  ($wavelength_acceptor - $wavelength_width/2, $wavelength_acceptor + $wavelength_width/2);
  my $filename = 0; #Raw spectrum filenames increase incrimentally and numerically starting at '0'.
                    #('0' filename is the spectra before any titrant has been added.)
  foreach (den($self)) {
    open_spectrum($self, $filename.chr(46).file_spectrum_extension($self, $ext), $delimiter);
    my @raw_spectrum = sort { $a->[0] <=> $b->[0] } @{ [raw_spectrum($self)] };
    my ($avg_sig, $avg_donor_sig, $avg_acceptor_sig, $avg_fret_sig) = (0, 0, 0, 0);
    my ($cnt, $total_sig) = (0, 0);
    my ($donor_cnt, $total_donor_sig) = (0, 0);
    my ($acceptor_cnt, $total_acceptor_sig) = (0, 0);
    for ( my $i = 0; $i < scalar(@raw_spectrum); $i++ ) {
      my ($wave, $sig) = ($raw_spectrum[$i][0], $raw_spectrum[$i][1]) ;
      #fluorescence_mode Either 'SINGLE_WAVELENGTH' or 'FRET'    
      if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') { 
        if (($wave >= $start) && ($wave <= $end)) {
	  $cnt++;
	  $total_sig += $sig;
        };
      }
      elsif ($fluorescence_mode eq 'FRET') {
        if (($wave >= $donor_start) && ($wave <= $donor_end)) {
	  $donor_cnt++;
	  $total_donor_sig += $sig;
        };
        if (($wave >= $acceptor_start) && ($wave <= $acceptor_end)) {
	  $acceptor_cnt++;
	  $total_acceptor_sig += $sig;
        };
      };  
    };
    if ($cnt) { $avg_sig = $total_sig/$cnt };
    if ($donor_cnt) { $avg_donor_sig = $total_donor_sig/$donor_cnt };
    if ($acceptor_cnt) { $avg_acceptor_sig = $total_acceptor_sig/$acceptor_cnt };
    if ($avg_acceptor_sig) { $avg_fret_sig =  $avg_donor_sig/$avg_acceptor_sig};
    if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') {
      push @signal, $avg_sig;
      print "Average signal at wavelength $wavelength nm is $avg_sig\n" if verbose($self);
     }
     elsif ($fluorescence_mode eq 'FRET') {
       push @signal, $avg_fret_sig;
       print "Average FRET signal at donor wavelength $wavelength_donor nm and acceptor wavelength $wavelength_acceptor nm is $avg_fret_sig\n" if verbose($self);
      };
      $filename++; #next integer increment filename 
  };
  signal($self, @signal);
  stitch_data($self);
  return $self
}

#Extract single wavelength fluorescence and FRET signals from one large file where the titration points
#are pairs of columns in one large CSV ASCII text file
#odd-numbered columns are wavelength and even-numbered columns are corresponding signals
sub extract_batch_wavelength_signal {
  my ($self, $file, $header, $fluorescence_mode, $delimiter) = @_;
  if ($fluorescence_mode) { fluorescence_mode($self, $fluorescence_mode) } else { $fluorescence_mode = fluorescence_mode($self) }
  $delimiter = chr(44) unless $delimiter; #default delimiter for batch is comma
  my ($spectra, $signal, $raw_data) = ([], [], []);
  print "fluorescence_mode is $fluorescence_mode\n"           if verbose($self);
  my $wavelength = wavelength($self);
  print "wavelength is $wavelength nm\n"                      if (verbose($self) && $wavelength);
  my ($wavelength_width, $wavelength_donor, $wavelength_acceptor) = (wavelength_width($self), wavelength_donor($self), wavelength_acceptor($self));
  print "wavelength_width is $wavelength_width nm\n"          if verbose($self);
  print "wavelength_donor is $wavelength_donor nm\n"          if (verbose($self) && $wavelength_donor);
  print "wavelength_acceptor is $wavelength_acceptor nm\n"    if (verbose($self) && $wavelength_acceptor);
  #start and end for averaging window are plus or minus half the width
  my ($start, $end) = ($wavelength - $wavelength_width/2, $wavelength + $wavelength_width/2);
  my ($donor_start, $donor_end) =  ($wavelength_donor - $wavelength_width/2, $wavelength_donor + $wavelength_width/2);
  my ($acceptor_start, $acceptor_end) =  ($wavelength_acceptor - $wavelength_width/2, $wavelength_acceptor + $wavelength_width/2);  
  open(FH, '<', $file) || die "Cannot open $file in extract_batch_wavelength_signal. $!";
  if ($header) { foreach (1..$header) { my $line = <FH> } }; #parse through header
  #open file contents into raw data, which is an array of arrays
  print "Opening $file for extract_batch_wavelength_signal.\n" if verbose($self); 
  my $end_of_data_boolean = 0;
  while (!$end_of_data_boolean) {
    my $line = <FH>; 
    $line = ts($line); #remove linefeed/carriage return whitespace from end of line
    $line =~ s/$delimiter$//; #remove free delimiter (often found as a comma in csv file) from end of line
    $end_of_data_boolean = 1 unless $line; #end of data signified with a blank line
    push @$raw_data, [ split(/$delimiter/, $line) ] if $line;
  };
  close(FH);
  #separate raw data into individual spectra
  my $number_columns = scalar(@{ $raw_data -> [0] });
  for (my $i = 0; $i < $number_columns; $i += 2) {
    my ($wavelength_column, $signal_column) = ($i, $i + 1);
    my $spectrum = [];
    foreach my $row (@$raw_data) { push @$spectrum, [ $row -> [$wavelength_column], $row -> [$signal_column] ] };
    push @$spectra, $spectrum;
  };
  #Now go through each spectrum in spectra array of arrays of arrays and get signals at single wavelength or FRET wavelengths
  foreach my $spectrum (@$spectra) {
    my ($avg_sig, $avg_donor_sig, $avg_acceptor_sig, $avg_fret_sig) = (0, 0, 0, 0);
    my ($cnt, $total_sig) = (0, 0);
    my ($donor_cnt, $total_donor_sig) = (0, 0);
    my ($acceptor_cnt, $total_acceptor_sig) = (0, 0);
    for ( my $i = 0; $i < scalar(@$spectrum); $i++ ) {
      my ($wave, $sig) = ($spectrum->[$i]->[0], $spectrum->[$i]->[1]) ;
      #fluorescence_mode Either 'SINGLE_WAVELENGTH' or 'FRET'    
      if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') { 
        if (($wave >= $start) && ($wave <= $end)) {
	  $cnt++;
	  $total_sig += $sig;
        };
      }
      elsif ($fluorescence_mode eq 'FRET') {
        if (($wave >= $donor_start) && ($wave <= $donor_end)) {
	  $donor_cnt++;
	  $total_donor_sig += $sig;
        };
        if (($wave >= $acceptor_start) && ($wave <= $acceptor_end)) {
	  $acceptor_cnt++;
	  $total_acceptor_sig += $sig;
        };
      };  
    };
    if ($cnt) { $avg_sig = $total_sig/$cnt };
    if ($donor_cnt) { $avg_donor_sig = $total_donor_sig/$donor_cnt };
    if ($acceptor_cnt) { $avg_acceptor_sig = $total_acceptor_sig/$acceptor_cnt };
    if ($avg_donor_sig) { $avg_fret_sig =  $avg_acceptor_sig/$avg_donor_sig}; 
    if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') {
      push @$signal, $avg_sig;
      print "Average signal at wavelength $wavelength nm is $avg_sig\n" if verbose($self);
     }
     elsif ($fluorescence_mode eq 'FRET') {
       push @$signal, $avg_fret_sig;
       print "Average FRET signal at donor wavelength $wavelength_donor nm and acceptor wavelength $wavelength_acceptor nm is $avg_fret_sig\n" if verbose($self);
      };
  };
  signal($self, @$signal);
  stitch_data($self);
  return $self
}

#Extract single wavelength fluorescence and FRET signals from one large file where the titration points
#are columns in one large CSV ASCII text file
#first column is the wavelength and all columns thereafter are the signals corresponding to each titration point
#Basically fluorescence intensity signals are the z-axis in a 3-D plot of the CSV data file
sub extract_z_batch_wavelength_signal {
  my ($self, $file, $header, $fluorescence_mode, $delimiter) = @_;
  if ($fluorescence_mode) { fluorescence_mode($self, $fluorescence_mode) } else { $fluorescence_mode = fluorescence_mode($self) }
  $delimiter = chr(44) unless $delimiter; #default delimiter for batch is comma
  my ($spectra, $signal, $raw_data) = ([], [], []);
  print "fluorescence_mode is $fluorescence_mode\n"           if verbose($self);
  my $wavelength = wavelength($self);
  print "wavelength is $wavelength nm\n"                      if (verbose($self) && $wavelength);
  my ($wavelength_width, $wavelength_donor, $wavelength_acceptor) = (wavelength_width($self), wavelength_donor($self), wavelength_acceptor($self));
  print "wavelength_width is $wavelength_width nm\n"          if verbose($self);
  print "wavelength_donor is $wavelength_donor nm\n"          if (verbose($self) && $wavelength_donor);
  print "wavelength_acceptor is $wavelength_acceptor nm\n"    if (verbose($self) && $wavelength_acceptor);
  #starts and ends for the signal averaging windows are plus or minus half the wavelength_width
  my ($start, $end) = ($wavelength - $wavelength_width/2, $wavelength + $wavelength_width/2);
  my ($donor_start, $donor_end) =  ($wavelength_donor - $wavelength_width/2, $wavelength_donor + $wavelength_width/2);
  my ($acceptor_start, $acceptor_end) =  ($wavelength_acceptor - $wavelength_width/2, $wavelength_acceptor + $wavelength_width/2);  
  open(FH, '<', $file) || die "Cannot open $file in extract_z_batch_wavelength_signal. $!";
  if ($header) { foreach (1..$header) { my $line = <FH> } }; #parse through header
  #open file contents into raw data, which is an array of arrays
  print "Opening $file for extract_z_batch_wavelength_signal.\n" if verbose($self); 
  my $end_of_data_boolean = 0;
  while (!$end_of_data_boolean) {
    my $line = <FH>; 
    $line = ts($line); #remove linefeed/carriage return whitespace from end of line
    $line =~ s/$delimiter$//; #remove free delimiter (comma often found in csv file) from end of line
    $end_of_data_boolean = 1 unless $line; #end of data signified with a blank line
    push @$raw_data, [ split(/$delimiter/, $line) ] if $line;
  };
  close(FH);
  #separate raw data into individual spectra
  #in the larger CSV file, the first column is the wavelength column and columns after it are the signal at each denaturant titration condition
  my $number_columns = scalar(@{ $raw_data -> [0] });
  for (my $i = 1; $i < $number_columns; $i++) {
    my ($wavelength_column, $signal_column) = (0, $i);
    my $spectrum = [];
    foreach my $row (@$raw_data) { push @$spectrum, [ $row -> [$wavelength_column], $row -> [$signal_column] ] };
    push @$spectra, $spectrum;
  };
  #Now go through each spectrum in spectra array of arrays of arrays and get signals at single wavelength or FRET wavelengths
  foreach my $spectrum (@$spectra) {
    my ($avg_sig, $avg_donor_sig, $avg_acceptor_sig, $avg_fret_sig) = (0, 0, 0, 0);
    my ($cnt, $total_sig) = (0, 0);
    my ($donor_cnt, $total_donor_sig) = (0, 0);
    my ($acceptor_cnt, $total_acceptor_sig) = (0, 0);
    for ( my $i = 0; $i < scalar(@$spectrum); $i++ ) {
      my ($wave, $sig) = ($spectrum -> [$i] -> [0], $spectrum -> [$i] -> [1]) ;
      #fluorescence_mode Either 'SINGLE_WAVELENGTH' or 'FRET'    
      if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') { 
        if (($wave >= $start) && ($wave <= $end)) {
	  $cnt++;
	  $total_sig += $sig;
        };
      }
      elsif ($fluorescence_mode eq 'FRET') {
        if (($wave >= $donor_start) && ($wave <= $donor_end)) {
	  $donor_cnt++;
	  $total_donor_sig += $sig;
        };
        if (($wave >= $acceptor_start) && ($wave <= $acceptor_end)) {
	  $acceptor_cnt++;
	  $total_acceptor_sig += $sig;
        };
      };  
    };
    if ($cnt) { $avg_sig = $total_sig/$cnt };
    if ($donor_cnt) { $avg_donor_sig = $total_donor_sig/$donor_cnt };
    if ($acceptor_cnt) { $avg_acceptor_sig = $total_acceptor_sig/$acceptor_cnt };
    if ($avg_donor_sig) { $avg_fret_sig =  $avg_acceptor_sig/$avg_donor_sig}; 
    if ($fluorescence_mode eq 'SINGLE_WAVELENGTH') {
      push @$signal, $avg_sig;
      print "Average signal at wavelength $wavelength nm is $avg_sig\n" if verbose($self);
     }
     elsif ($fluorescence_mode eq 'FRET') {
       push @$signal, $avg_fret_sig;
       print "Average FRET signal at donor wavelength $wavelength_donor nm and acceptor wavelength $wavelength_acceptor nm is $avg_fret_sig\n" if verbose($self);
      };
  };
  signal($self, @$signal);
  stitch_data($self);
  return $self
}

#Stitch single-wide denaturant or ligand and signal arrays into a two-wide array of arrays called 'data',
#where column one (x-axis) is denaturant or ligand concentration and column two (y-axis) is corresponding signal.
# the module property titration_mode decides whether column x is denaturant or ligand
sub stitch_data {
  my $self = shift;
  my @data = ();
  if (titration_mode($self) eq 'DENATURANT') {
    for (0..(scalar(den($self)) - 1)) { push @data, [ [den($self)]->[$_], [signal($self)]->[$_] ] };
  }  
  elsif (titration_mode($self) eq 'LIGAND') {
    for (0..(scalar(ligand($self)) - 1)) { push @data, [ [ligand($self)]->[$_], [signal($self)]->[$_] ] };
  };
  data($self, @data); #store data in module 
  return $self
}

#Saves the 2-column data with denaturant or ligand concentration in the x-axis column and
# spectroscopic signal in the y-axis column
# the module property titration_mode decides whether column x is denaturant or ligand
sub save_data {
  my ($self, $out_file, $delimiter, $data_array_ref) = @_;
  $out_file  = file_out($self) unless $out_file;
  $delimiter = delimiter($self, chr(9)) unless $delimiter;
  delimiter($self, $delimiter);
  my @data = ();
  if (ref($data_array_ref) eq 'ARRAY') { @data = @$data_array_ref } else { @data = data($self) }; 
  open(FH, '>', file_out($self, $out_file)) || die "Cannot open $out_file in save_data. $!"; 
  my $cr = chr(13).chr(10); #carriage return and linefeed for PC and Mac compatibility 
  if (header($self) && !$data_array_ref) {
    if (titration_mode($self) eq 'DENATURANT') {
      print FH den_name($self).$delimiter.signal_name($self).$cr;
    }
    elsif (titration_mode($self) eq 'LIGAND') {
      print FH ligand_name($self).$delimiter.signal_name($self).$cr
    };
  };
  if (titration_mode($self) eq 'DENATURANT') {
    print den_name($self).$delimiter.signal_name($self).$cr if (verbose($self) && header($self) && !$data_array_ref); 
  }
  elsif (titration_mode($self) eq 'LIGAND') {
    print ligand_name($self).$delimiter.signal_name($self).$cr if (verbose($self) && header($self) && !$data_array_ref); 
  };
  foreach my $datum (@data) {
    print FH $datum->[0].$delimiter.$datum->[1].$cr;
    print $datum->[0].$delimiter.$datum->[1].$cr if verbose($self);
  };
  close(FH);
  return $self
}

#trim spaces - remove trailing white spaces, i.e. the carriage return and linefeed for PC/DOS compatible file data
sub ts { my $line = shift; $line =~ s/\s+$//; return $line }

############################# PROPERTY SETS/GETS ###################################

### ARRAY PROPERTIES

#Array of arrays used in timed titrations where the time (column zero) and signal (column one) are given
sub raw_data {
 my $self = shift;
 if (@_) { @{ $self->{RAW_DATA} } = @_ };
 return @{ $self->{RAW_DATA} };
}

#The currently opened spectrum file as an array of arrays given in wavelength (column zero) and signal (column one) format
sub raw_spectrum {
  my $self = shift;
  if (@_) { @{ $self->{RAW_SPECTRUM} } = @_ }
  return @{ $self->{RAW_SPECTRUM} };
}

#Array containing the list of titrant volumes per point
sub titration {
  my $self = shift;
  if (@_) { @{ $self->{TITRATION} } = @_ };
  return @{ $self->{TITRATION} };
}

#denaturant concentration per step including the starting cuvette concentration
sub den {
  my $self = shift;
  if (@_) { @{ $self->{DEN} } = @_ };
  return @{ $self->{DEN} };
}

#ligand concentrations per step inclusing the starting cuvette concentration
sub ligand {
  my $self = shift;
  if (@_) { @{ $self->{LIGAND} } = @_ };
  return @{ $self->{LIGAND} };
}

#Array of arrays of the volumes delivered per ligand stock concentration 
sub ligand_titration {
  my $self = shift;
  if (@_) { @{ $self->{LIGAND_TITRATION} } = @_ };
  return @{ $self->{LIGAND_TITRATION} };
}

#Array of the ligand stock concentrations in the order used in a ligand titration
sub ligand_stock_concentrations {
  my $self = shift;
  if (@_) { @{ $self->{LIGAND_STOCK_CONCENTRATIONS} } = @_ };
  return @{ $self->{LIGAND_STOCK_CONCENTRATIONS} };
}   

#corresponding signal array for each titration step including the signal at the starting cuvette concentration
sub signal {
  my $self = shift;
  if (@_) { @{ $self->{SIGNAL} } = @_ };
  return @{ $self->{SIGNAL} };
}

#Final post-processing array of arrays of the denaturant or ligand concentration in column zero
#and the corresponding spectroscopic signal in column one
sub data {
  my $self = shift;
  if (@_) { @{ $self->{DATA} } = @_ }; 
  return @{ $self->{DATA} };
}

### SCALAR PROPERTIES

#TITRATION_PROGRAM sets the titration protocol used, for example, two segments are as follows
# in its shorthand form, i.e., 30_25_48_64.
# The first segment has 30 pts of 25 uL titrant volume; the second segment has 48 pts of 64 uL titrant volume.

# In a table it would be:
#
# SEG. 	PTS.	VOLUME (UL)
# ============================
# SEG1	30	25
# SEG2	48	64
# ============================
#
# NOTE the TITRATION_PROGRAM must separate the values for the segments using only an underscore
# and no spaces or other characters—only numbers and underscores.
# A simple one-segment program would be 153_25, which is 153 points of 25 uL titrant volume per point
sub titration_program {
  my $self = shift;
  if (@_) {$self->{TITRATION_PROGRAM} = shift };
  return $self->{TITRATION_PROGRAM};
}

#volume of the cuvette in uL
sub cuvette_volume {
  my $self = shift;
  if (@_) { $self->{CUVETTE_VOLUME} = shift };
  return $self->{CUVETTE_VOLUME};
}

#conentration in molar of denaturant in the titrant in molar
sub titrant_conc {
  my $self = shift;
  if (@_) { $self->{TITRANT_CONC} = shift };
  return $self->{TITRANT_CONC};
}

#initial denaturant concentration in molar of the cuvette
sub cuvette_conc {    
  my $self = shift;
  if (@_) { $self->{CUVETTE_CONC} = shift };
  return $self->{CUVETTE_CONC};
}

#initial ligand concentration in the cuvette
sub cuvette_ligand_concentration {
  my $self = shift;
  if (@_) { $self->{CUVETTE_LIGAND_CONCENTRATION} = shift };
  return $self->{CUVETTE_LIGAND_CONCENTRATION};
}

#Either 'DENATURANT' or 'LIGAND' titration experiment
sub titration_mode {
  my $self = shift;
  if (@_) { $self->{TITRATION_MODE} = shift };
  return $self->{TITRATION_MODE};
}

#Either 'SINGLE_WAVELENGTH' or 'FRET' experiment
sub fluorescence_mode {
  my $self = shift;
  if (@_) { $self->{FLUORESCENCE_MODE} = shift};
  return $self->{FLUORESCENCE_MODE};
}

#When conjoining two timed titration files, this is the time offset added to the entire second file
sub time_offset {
  my $self = shift;
  if (@_){ $self->{TIME_OFFSET} = shift };
  return $self->{TIME_OFFSET};
}

#length of time at the beginning a titration before any titrant has been delivered
sub time_init {
  my $self = shift;
  if (@_) { $self->{TIME_INIT} = shift };
  return $self->{TIME_INIT};
}

#length of of time per titration step in a timed titration experiment
sub time_per_bin {
  my $self = shift;
  if (@_) { $self->{TIME_PER_BIN} = shift };
  return $self->{TIME_PER_BIN};
}

#Used to mark start of timed denaturant titration experiment
sub time_start {
  my $self = shift;
  if (@_) { $self->{TIME_START} = shift };
  return $self->{TIME_START};
}

#End of time in timed denaturant titration experiment. Used to mark end of experiment
sub time_end {
  my $self = shift;
  if (@_) { $self->{TIME_END} = shift };
  return $self->{TIME_END};
}

#wavelength is the wavelength used in signal extraction using SINGLE_WAVELENGTH fluorescence mode
sub wavelength {
  my $self = shift;
  if (@_) { $self->{WAVELENGTH} = shift };
  return $self->{WAVELENGTH};
}

#like slit width or averaging window for measuring signal at wavelength, wavelength_donor, and wavelength_acceptor  
sub wavelength_width {
  my $self = shift; 
  if (@_) { $self->{WAVELENGTH_WIDTH} = shift };
  return $self->{WAVELENGTH_WIDTH};
}

#wavelength in nm of the donor in a FRET experiment
sub wavelength_donor {
  my $self = shift;
  if (@_) { $self->{WAVELENGTH_DONOR} = shift };
  return $self->{WAVELENGTH_DONOR};
}

#wavelength in nm of the acceptor in a FRET experiment
sub wavelength_acceptor {
  my $self = shift;
  if (@_) { $self->{WAVELENGTH_ACCEPTOR} = shift };
  return $self->{WAVELENGTH_ACCEPTOR};
}

#filename of timed titration raw data file
sub file {
  my $self = shift;
  if (@_) { $self->{FILE} = shift };
  return $self->{FILE};
}

#filename of currently opened spectra file
sub file_spectrum {
  my $self = shift;
  if (@_) { $self->{FILE_SPECTRUM} = shift };
  return $self->{FILE_SPECTRUM};
}

#file extention used in opening a batch of spectra files
sub file_spectrum_extension {
  my $self = shift;
  if (@_) { $self->{FILE_SPECTRUM_EXTENSION} = shift };
  return $self->{FILE_SPECTRUM_EXTENSION};
}

#length in number of lines of header information in spectra files
#the program will skip over the header info when loading a spectrum
sub file_spectrum_header_length {
  my $self = shift;
  if (@_) { $self->{FILE_SPECTRUM_HEADER_LENGTH} = shift }
  return $self->{FILE_SPECTRUM_HEADER_LENGTH};
}

#input parameter filename with hash keys to the Den.pm module and their respective values
sub parameter_file {
  my $self = shift;
  if (@_) { $self->{PARAMETER_FILE} = shift };
  return $self->{PARAMETER_FILE};
}

#output filename
sub file_out {
  my $self = shift;
  if (@_) { $self->{FILE_OUT} = shift };
  return $self->{FILE_OUT};
}

#file input and output data separator
sub delimiter {
  my $self = shift;
  if (@_) { $self->{DELIMITER} = shift };
  return $self->{DELIMITER};
}

#name of denaturant used, e.g. 'Urea' or 'Gdm'
sub den_name {
  my $self = shift;
  if (@_) { $self->{DEN_NAME} = shift };    
  return $self->{DEN_NAME};
}

#Name of the ligand used in a ligand titration
sub ligand_name {
  my $self = shift;
  if (@_) { $self->{LIGAND_NAME} = shift };    
  return $self->{LIGAND_NAME};
}

#Name of signal used, e.g. 'FL', 'CD', or 'FRET'
sub signal_name {
  my $self = shift;
  if (@_) { $self->{SIGNAL_NAME} = shift };
  return $self->{SIGNAL_NAME};
}

#header is a Perl boolean. If true (set to non-zero value or string) adds header of denaturant name and signal name
#to output data file
#False value is zero or undef
sub header {
  my $self = shift;
  if (@_) { $self->{HEADER} = shift };
  return $self->{HEADER};
}

#verbose is a Perl boolean. If true (set to non-zero value or string), then
#module functions print their processing info to the cmd prompt
#False value is zero or undef
sub verbose {
  my $self = shift;
  if (@_) { $self->{VERBOSE} = shift };
  return $self->{VERBOSE};
}

1; #Must return 1 at end of module
