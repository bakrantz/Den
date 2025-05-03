#####################################################################################
#################                   Hoh.pm                           ################
#################             Krantz Lab       July 2008             ################
#####################################################################################
#################  Module to manipulate datasets in hash format      ################
#####################################################################################
package Hoh;
################################# THE CONSTRUCTOR ###################################
sub new {
  my $proto                            = shift;
  my $package                          = ref($proto) || $proto;
  my $self                             = {};                      #the anon hash
     $self->{HOH}                      = {};
     $self->{SEL}                      = {};
     $self->{PRINT_ORDER}              = [];
     $self->{PRINT_ORDER_FILE}         = undef;
     $self->{FILE_OUT}                 = undef;
     $self->{FILE}                     = undef;  
     $self->{KEY_NAMES}                = undef; 
     $self->{DELIMITER}                = undef; 
     $case->{CASE_SENSITIVE}           = undef;
     $self->{SORT_BY}                  = undef;
     $self->{HEADER_OFF}               = undef;
     $self->{NULL}                     = '--';
  return bless($self, $package);                                  #return thy self
}
############################### METHOD SUBS ########################################
sub join_files {
 my ($self, $file_a, $file_b, $file_out, $file_print_order) = (@_);
 join_hohs($self, +{ file_to_hoh($self, $file_a) }, +{ file_to_hoh($self, $file_b) });
 $file_out = file_out($self) unless $file_out;
 $file_out = file_out($self, 'outfile.txt') unless file_out($self);
 if ($file_print_order) { open_print_order_file($self, $file_print_order) }
 elsif (print_order_file($self)) { open_print_order_file($self, print_order_file($self)) };
 unless (print_order($self)) { print_order($self, ( sort { $a cmp $b } keys %{ +{ hoh($self) }->{ shift(@{ [keys %{ +{ hoh($self) } }] }) } } ) ) };
 save($self, file_out($self, $file_out), delimiter($self), sort_by($self), header_off($self));
 return $self;
};

sub save {
 my ($self, $file_out, $d, $sort_by, $header_off) = @_;
 $file_out = file_out($self) unless $file_out;
 $d = delimiter($self) unless $d;
 $sort_by = sort_by($self) unless $sort_by;
 $header_off = header_off($self) unless $header_off;
 scalar_to_file( hoh_to_scalar($self, +{ hoh($self) }, [print_order($self)], sort_by($self, $sort_by), abs(header_off($self, $header_off) - 1), delimiter($self, $d)), file_out($self, $file_out));
 return $self;
}

sub open_print_order_file { print_order($_[0], (file_to_array( print_order_file($_[0], $_[1]) ) ) ); return $_[0]; }

sub join_hohs {
 my ($self, $i, $j) = @_;
 foreach my $r (rows($self)) { foreach my $c (cols($self)) { $i->{$r}->{$c} = $j->{$r}->{$c} } };
 hoh($self, %$i);
 return $self;
}

sub save_as_chimera_attribute {
 my ($self, $attr_name, $col_attr, $col_resi) = (@_);
 $attr_name =~ s/[\W_]+//g;
 $attr_name = 'attr'.$attr_name if ($attr_name =~ /^\d/); 
 $attr_name = lc(substr($attr_name, 0, 1)).substr($attr_name, 1, length($attr_name) - 1); 
 my $cr = chr(13).chr(10);
 print "Making a chimera attribute file from a Hoh text file...\nHoh file: ".file($self)."\n".
       "Chimera attribute name: $attr_name\nColumn in Hoh file to make attribute from: $col_attr\n";
 open(FH, ">$attr_name.txt");
 print FH join($cr, @ {[ "# From Krantz Lab Hoh.pm", "# Chimera attribute output file for residues", 
  "# Text file is ".file($self), "attribute: $attr_name", "match mode: 1-to-1", "recipient: residues" ] }).$cr;
 $col_resi = undef if ($col_resi eq key_names($self)); 
 $col_resi = undef unless ( exists +{ hoh($self) }->{ shift( @{ [rows($self)] } ) }->{$col_resi} );
 if ($col_resi) { foreach (rows($self)) { print FH "\t:".+{ hoh($self) }->{$_}->{$col_resi}."\t".+{ hoh($self) }->{$_}->{$col_attr}.$cr } }
 else { foreach (rows($self)) { print FH "\t:$_\t".+{ hoh($self) }->{$_}->{$col_attr}.$cr } };
 close(FH);
 print "Finished writing $attr_name.txt\n";
 return $self;
 }

sub file_to_array { 
 my @a = split(/\s+/, file_to_scalar($_[0]));
 for my $i (0..$#a) { $a[$i] = cs($self, $a[$i]) };
 return @a;
}

sub scalar_to_file { open (FH, ">$_[1]") || die "Cannot open file. Aborting. $!"; print FH "$_[0]"; }

sub file_to_scalar {
 open(FH, "<$_[0]") || die "$0 can't open $_[0] using file_to_scalar method. $!\n";
 return do { local $/;  <FH> }; #/change the default file separator locally 
}

sub hoh_to_scalar {
 my ($self, $hoh, $print_order, $col, $head, $d) = (@_);
 my ($cr, $scaler, $print) = (chr(13).chr(10), undef, []);
 $d = delimiter($self) unless $d;
 $d = delimiter($self, "\t") unless $d;
 my @rows  = sort { $a cmp $b } keys %$hoh;
   print "The print order in hoh_to_scalar is...\n".join("\n", @$print_order);
 @$print_order = sort { $a cmp $b } keys %$href unless @$print_order;
 foreach my $it (@$print_order) { push @$print, $it if exists $hoh->{$rows[0]}->{$it} };
 $scalar = key_names($self).$d.join("$d", @$print).$cr if $head; #head: boolean for header printing
 @rows = sort { $hoh->{$a}->{$col} <=> $hoh->{$b}->{$col} } keys %$hoh if ($hoh->{$rows[0]}->{$col});
 foreach my $r (@rows) {
  $scalar .= $r.$d;
  foreach (@$print) { $scalar .= $hoh->{$r}->{$_}.$d };
  $scalar =~ s/$d$/$cr/;
 };
 return $scalar;
};

sub open_hoh { my $self = shift; file_to_hoh($self, @_); return $self; }

sub file_to_hoh {
 my ($self, $file, $d, $cs) = @_; #Note that the files are in a key-tab-value format
 my ($tab, $hoh)    = (chr(9), {}); 
 $file = file($self) unless $file;
 $d    = delimiter($self) unless $d;
 $cs   = case_sensitive($self) unless $cs; 
 open(FH, "<", file($self, $file)) || die "$0 can't open $file using file_to_hoh method. $!";
 my $line =  <FH>;
 my @keys = line_to_array(cs($self, ts($line), $cs), delimiter($self, $d));
 while ($line = <FH>) {
  my @data = line_to_array(cs($self, ts($line)), $d);
  for (1..$#data) { $hoh -> {$data[0]} -> {$keys[$_]} = $data[$_] };
 };
 close(FH);
 key_names($self, $keys[0]);
 return hoh($self, %$hoh);
}

sub copy {
 my ($self, $cols, $rows) = @_;
 my %hoh = hoh($self);
 hoh($self, col($self, $cols));
 my %sel = row($self, $rows);
 hoh($self, %hoh);
 return sel($self, %sel);
}

sub cut {
 my ($self, $cols, $rows) = @_;
 my %cut = copy($self, $cols, $rows);
 my %hoh = hoh($self);
 hoh($self, %cut); 
 foreach my $rr (rows($self)) { foreach my $cc (cols($self)) { clear($self, $cc, $rr) } };
 join_hohs($self, \%hoh, +{ hoh($self) });
 return sel($self, %cut);
}

sub paste { return join_hohs($_[0], +{hoh($self)}, +{sel($self)}) }

sub clear_clipboard {
 my $self = shift;
 my %hoh = hoh($self);
 hoh($self, sel($self)); 
 foreach my $rr (rows($self)) { foreach my $cc (cols($self)) { clear($self, $cc, $rr) } };
 hoh($self, %hoh);
 return $self;
}

sub save_sel { save_clipboard(@_) }

sub save_clipboard {
 my ($self, $file_out, $d, $sort_by, $header_off) = @_;
 $file_out = file_out($self) unless $file_out;
 $d = delimiter($self) unless $d;
 $sort_by = sort_by($self) unless $sort_by;
 $header_off = header_off($self) unless $header_off;
 scalar_to_file( hoh_to_scalar($self, +{ sel($self) }, [], sort_by($self, $sort_by), abs(header_off($self, $header_off) - 1), delimiter($self, $d)), file_out($self, $file_out));
 return $self;
}

sub del {
 my ($self, $key) = @_;
 my $hoh = +{ hoh($self) };
 if    (is($key, [cols($self)])) { foreach (rows($self)) { delete $hoh->{$_}->{$key} } }
 elsif (is($key, [rows($self)])) { delete $hoh->{$key} };
 hoh($self, %$hoh);
 return $self;
}

sub val {
 my ($self, $c, $r, $val) = @_;
 my ($cc, $rr) = (undef, undef);
 my $hoh = +{ hoh($self) };
 $cc = $c if is($c, [cols($self)]);
 $rr = $r if is($r, [rows($self)]);
 if ($val) { $hoh->{$rr}->{$cc} = $val if ($cc && $rr) }
 else      { $val = $hoh->{$rr}->{$cc} if ($cc && $rr) };
 hoh($self, %$hoh);
 return $val;
}

sub nullify { clear(@_) }

sub clear { val($_[0], $_[1], $_[2], null($self)) }

sub col {
 my $self = shift();
 my @c    = @_;
 my ($col, $columns) = ([], {});
 if    (scalar(@c) == 1) {
  if (ref(\$c[0]) eq 'SCALAR') { push @$col, $c[0] if is($c[0], [cols($self)]) };
  if (ref( $c[0]) eq 'ARRAY' ) { foreach (@{$c[0]}) { push @$col, $_ if is($_, [cols($self)]) } };
 }  
 elsif (scalar(@c) >  1) { for my $i (0..$#c) { push @$col, $c[$i] if is($c[$i], [cols($self)]) } }
 else                    { return undef };
 foreach my $cc (@$col)  { foreach my $rr (rows($self)) { $columns->{$rr}->{$cc} = +{ hoh($self) }->{$rr}->{$cc} } };
 return sel($self, %$columns);
}

sub row {
 my $self = shift();
 my @r    = @_;
 my ($row, $rows) = ([], {});
 if    (scalar(@c) == 1) {
  if (ref(\$r[0]) eq 'SCALAR') { push @$row, $r[0] if is($r[0], [rows($self)]) };
  if (ref( $r[0]) eq 'ARRAY' ) { foreach (@{$r[0]}) { push @$row, $_ if is($_, [rows($self)]) } };
 }  
 elsif (scalar(@r) >  1) { for my $i (0..$#r) { push @$row, $r[$i] if is($r[$i], [rows($self)]) } }
 else                    { return undef };
 foreach my $cc (@$col)  { foreach my $rr (rows($self)) { $columns->{$rr}->{$cc} = +{ hoh($self) }->{$rr}->{$cc} } };
 return sel($self, %$rows);
}

sub ts { my $in = shift; $in =~ s/\s+$//; return $in; }

sub cs {
 my ($self, $it, $cs) = @_;
 $cs = case_sensitive($self) unless $cs;
 $it = uc($it) unless case_sensitive($self, $cs);
 return $it;
}

sub line_to_array {
 if ($_[1]) { return split(/\s*$_[1]\s*/, $_[0]) if $_[0] =~ /\s*$_[1]\s*/ }
 else       { return split(/\s+/, $_[0]) if $_[0] =~ /\s+/ };
}

sub is {
 my ($is_something, $in, $here) = ($_[0], $_[1], {});
 if    (ref($in)  eq 'ARRAY' ) { foreach (@$in) { $here -> {$_} = 1 } }
 elsif (ref($in)  eq 'HASH'  ) { %$here = %$in                        }
 elsif (ref($in)  eq 'SCALAR') { $here -> {$$in} = 1                  }
 elsif (ref(\$in) eq 'SCALAR') { $here -> {$in}  = 1                  }
 else { return undef };
 return exists $here -> {$is_something}; 
}

############################# PROPERTY SETS/GETS ###################################

###ALL HASHES

sub hoh {
  my $self = shift;
  if (@_) { %{ $self->{HOH} } = @_ };
  return %{ $self->{HOH} };
}

sub sel {
  my $self = shift;
  if (@_) { %{ $self->{SEL} } = @_ };
  return %{ $self->{SEL} };
}

sub clipboard { sel(@_) }

sub print_order {
 my $self = shift;
 if (@_) {
  my @in =  @_;
  for my $i (0..$#in) { $in[$i] = cs($self, $in[$i]) };
  @{ $self->{PRINT_ORDER} } = @in; 
};
  return @{ $self->{PRINT_ORDER} };
}

### ALL ARRAY PROPERTIES ###

sub cols { return @{ [ sort { $a cmp $b } keys %{ +{ hoh($_[0]) } -> { shift( @{ [ keys %{ +{ hoh($_[0]) }  } ] } ) } } ] } }

sub rows { return @{ [ sort { $a cmp $b } keys %{ +{ hoh($_[0]) } } ] } }

###ALL SCALARS
sub ncols { return scalar(cols($_[0])) }

sub nrows { return scalar(rows($_[0])) }

sub sort_by {
 my $self = shift;
 if (@_) { $self->{SORT_BY} = shift };
 return $self->{SORT_BY};
}

sub print_order_file {
 my $self = shift;
 if (@_) { $self->{PRINT_ORDER_FILE} = shift }; 
 return $self->{PRINT_ORDER_FILE};
}

sub file_out {
 my $self = shift;
 if (@_) { $self->{FILE_OUT} = shift }; 
 return $self->{FILE_OUT};
}

sub file {
 my $self = shift;
 if (@_) { $self->{FILE} = shift }; 
 return $self->{FILE};
}

sub key_names {
 my $self = shift;
 if (@_) { $self->{KEY_NAMES} = shift }; 
 return $self->{KEY_NAMES};
}

sub delimiter {
 my $self = shift;
 if (@_) {  $self->{DELIMITER} = shift };
 return $self->{DELIMITER};
}

sub null { return $self->{NULL} }

sub case_sensitive {
 my $self = shift;
 if (@_) { $self->{CASE_SENSITIVE} = shift };
 return $self->{CASE_SENSITIVE};
}

sub header_off {
 my $self = shift;
 if (@_) { $self->{HEADER_OFF} = shift };
 return $self->{HEADER_OFF};
}

return 1;