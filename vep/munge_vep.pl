#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# input file is 5 columns
# 1. gene_id
# 2. population_id
# 3. consequence
# 4. homo/het
# 5. individuals (comma sep)

# generate a tab delimited file
# 1. id (aka gene_id)
# 2. capabilities (VEP)
# 3-n. VEP_<consequence>_<homo/het>_<population_id>_attr_ss (these can be combined at query time)
# n+1 VEP_merged_EMS_attr_ss
# n+1 VEP_merged_NAT_attr_ss

my @poptype = qw(NA EMS EMS NAT EMS NAT); # HARDCODED classifications of cabot:sorghum_bicolor_variation_6_87_30/population

my %hsh;
my %merged;
my %out;
while (<>) {
    chomp;
    my ($id, $pop_id, $conseq, $hh, $inds) = split /\t/, $_;
    $out{$id}{id} = $id;
    $out{$id}{capabilities} = "VEP";
    $hsh{$conseq}{$hh}{$pop_id}{$id} = $inds;
    for my $ind (split /,/, $inds) {
        $merged{$poptype[$pop_id]}{$id}{$ind}=1;
    }
}

my @fields = ('id','capabilities');
for my $c (keys %hsh) {
    for my $hh (keys %{$hsh{$c}}) {
        for my $p (keys %{$hsh{$c}{$hh}}) {
            my $field = join ("__", ("VEP",$c,$hh,$p,"attr_ss"));
            push @fields, $field;
            for my $id (keys %{$hsh{$c}{$hh}{$p}}) {
                $out{$id}{$field} = $hsh{$c}{$hh}{$p}{$id};
            }
        }
    }
}
for my $p (keys %merged) {
    my $field = "VEP__merged__${p}__attr_ss";
    push @fields, $field;
    for my $id (keys %{$merged{$p}}) {
        $out{$id}{$field} = join(',', sort keys %{$merged{$p}{$id}})
    }
}

print join("\t", @fields),"\n";
for my $id (keys %out) {
    print join("\t", map { $out{$id}{$_} || '' } @fields),"\n";
}
