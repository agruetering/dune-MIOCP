#! /usr/bin/perl

foreach $file (@ARGV) {
  open(FILE, $file);

  print "\% $file\n";
  if($file =~ /cut_(\d+)\.txt/) {
      		$nt = $1;
      }
  while (<FILE>) {
      if($_ =~ /^nodes:(\d+)/) {
	  $nodes = $1;
      }
      if($_ =~ /^dt\_(\d+):(\d+\.?\d*)/) {
	  $stepwidth[$1] = $2;
      }
      if($_ =~ /^u\_(\d+):(\d+\.?\d*)/) {
	  if($1 < 0 || $1 >= $nodes) {
	      print "error: index $1$ out of range {0,...,";
	      print $nodes[$nt]-1;
	      print "}\n";
	      exit(1);
	  }
	  $u[$1] = $2;
      }
  }
  close(FILE);

  $sum=0;
  for($k=0; $k<$nodes; $k++){
	 $sum+=$stepwidth[$k]; 
  }
  $endtime=$sum;
  
  my @cmd=('../build-cmake/src/robj');
  push @cmd, $nt;
  for($i=0; $i< $nodes; $i++){
	push @cmd, $u[$i];
  }
  system(@cmd);

}

