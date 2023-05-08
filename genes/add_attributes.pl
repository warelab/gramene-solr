#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use JSON;

my $attr_file = shift @ARGV;
-f $attr_file or die "usage: $0 attributes.txt < solr_genes.json";

# read attributes into memory
my %lut; # $lut{geneID}{attribute} = value;
# assume first line is a header with attribute names
open(my $fh, "<", $attr_file);
my $header = <$fh>;
chomp $header;
my ($id_col, @attrs) = split /\t/, $header;
my %is_multi;
for my $attr (@attrs) {
    $attr =~ m/^\w+_attr_[sif]s?$/ or $attr eq "capabilities" or die "error parsing $attr\n";
    $is_multi{$attr} = ($attr =~ m/_[sif]s$/);
}
while (<$fh>) {
    my ($id, @etc) = split /\t/, $_;
    scalar @etc == scalar @attrs or die "number of attributes (".scalar(@attrs).") differs in line (".scalar(@etc).")\n$_\n";
    chomp $etc[-1];
    for (my $i=0; $i<@etc; $i++) {
        my $attr = $attrs[$i];
        my $value = $etc[$i];
        next if ($value eq ''); 
        if ($is_multi{$attr}) {
            my @x = split /,/, $value;
            $value = \@x;
        }
        $lut{$id}{$attr} = $value;
    }
}
close $fh;

print STDERR "finished reading attributes into memory\n";

# process the json on stdin
# assume it is an array of objects with one object per line
while (<>) {
    if (/^{/) {
        my $doc = decode_json $_;
        my $id = $doc->{"id"};
        if ($lut{$id}) {
            for my $attr (keys %{$lut{$id}}) {
                if (not exists $doc->{$attr}) {
                    $doc->{$attr} = $lut{$id}{$attr};
                }
                else {
                    if ($attr eq "capabilities") {
                        push @{$doc->{$attr}}, $lut{$id}{$attr};
                    }
                }
            }
        }
        my $json = encode_json $doc;
        print "$json\n";
    }
    else {
        print $_;
    }
}

