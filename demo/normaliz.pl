#!/usr/bin/perl -w

open FILE,"rev1.bedGraph";
open OUT,">rev1.norm.bedGraph";

while(<FILE>){
	chomp;
	my ($chr,$start,$end,$val) = split;
	$val = sprintf("%.5f",$val/259);
	print OUT "$chr\t$start\t$end\t$val\n";

}
