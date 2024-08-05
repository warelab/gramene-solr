#!/usr/local/bin/perl 

=head1 NAME

document this code

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use autodie;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

get_consequences.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to dump
    --debug
    --nowrite don't actually update the database
    --trace output of TRaCE

=head1 OPTIONS

=over 4

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose transcripts are to be dumped

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=back

=head1 ARGUMENTS


=cut

my ($species, $registry);
my ($debug);

{ #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions(
    "help|?"=>\$help
    ,"man"=>\$man
	,"species=s"=>\$species
	,"registry=s"=>\$registry
	,"debug"=>\$debug
  ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
}


my %consequence_level = (
    "3_prime_UTR_variant"              => 0,
    "5_prime_UTR_variant"              => 0,
    coding_sequence_variant            => 0,
    frameshift_variant                 => 0,
    downstream_gene_variant            => 0,
    inframe_deletion                   => 0,
    inframe_insertion                  => 0,
    intron_variant                     => 0,
    missense_variant                   => 0,
    non_coding_transcript_exon_variant => 0,
    non_coding_transcript_variant      => 0,
    mature_miRNA_variant               => 0,
    protein_altering_variant           => 0,
    splice_acceptor_variant            => 2,
    splice_donor_variant               => 2,
    splice_region_variant              => 0,
    start_lost                         => 2,
    stop_gained                        => 2,
    stop_lost                          => 1,
    stop_retained_variant              => 0,
    synonymous_variant                 => 0,
    transcript_ablation                => 0,
    upstream_gene_variant              => 0
);

my (@these_cs, @not_these);
for my $c (keys %consequence_level) {
    if ($consequence_level{$c} > 0) {
        push @these_cs, $c;
    } else {
        push @not_these, $c;
    }
}
my $in_these_cs = join (",", map {"'$_'"} @these_cs);
my $not_these_cs = join (",", map {"'$_'"} @not_these);

# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );

# connect to core db
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );
my $gdbc = $ENS_DBA->dbc->db_handle;

# connect to variation DB
my $ENS_VAR_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'variation' );
$ENS_VAR_DBA || pod2usage( "\nNo variation DB for $species set in $registry\n" );
my $dbc = $ENS_VAR_DBA->dbc->db_handle;

# get canonical transcript to gene lookup table
my %t2g;
my $sth = $gdbc->prepare(qq{
    SELECT t.stable_id, g.stable_id
    FROM transcript t, gene g
    WHERE t.transcript_id = g.canonical_transcript_id and g.is_current = 1
});
$sth->execute();
my ($tr_id, $g_id);
$sth->bind_columns(\$tr_id,\$g_id);
$t2g{$tr_id} = $g_id while $sth->fetch();
$sth->finish();

# get the taxonomy id
$sth = $gdbc->prepare("select meta_value from meta where meta_key = 'species.production_name' and species_id = 1");
$sth->execute();
my ($prod_name);
$sth->bind_columns(\$prod_name);
$sth->fetch();
$sth->finish();
print STDERR "prod_name $prod_name\n";

# get seq_regions
$sth = $dbc->prepare(qq{
  SELECT DISTINCT(seq_region_id)
  FROM variation_feature
});
$sth->execute();

my ($sr_id, @sr_ids);
$sth->bind_columns(\$sr_id);
push @sr_ids, $sr_id while $sth->fetch();
$sth->finish();

# get a lookup table for genotypes
$sth = $dbc->prepare(qq{
  SELECT gc.genotype_code_id, ac.allele, gc.haplotype_id from genotype_code gc, allele_code ac where gc.allele_code_id = ac.allele_code_id
});
$sth->execute();
my %genotype_to_ACGT;
while (my ($gcode,$acgt,$hap) = $sth->fetchrow_array) {
  $genotype_to_ACGT{$gcode}{$acgt}=$hap;
}
$sth->finish();

# get a lookup table mapping individual_id to name, source and population
$sth = $dbc->prepare(qq{
  SELECT s.sample_id, s.name, st.source_id, sp.population_id, st.study_type
  from sample s, study st, sample_population sp
  where s.sample_id = sp.sample_id and s.study_id = st.study_id
});
# $sth = $dbc->prepare(qq{
#   SELECT s.sample_id, s.individual_id, s.name, st.source_id, st.study_type
#   from sample s, study st
#   where s.study_id = st.study_id
# });
$sth->execute();
my %individual_pop;
my %individual_source;
my %individual_name;
my %study_type;
while (my ($sample,$name,$source,$pop,$stype) = $sth->fetchrow_array) {
  $individual_source{$sample} = $source;
  $individual_name{$sample} = $name;
  $individual_pop{$sample} = $pop;
  $study_type{$pop} = ($stype and $stype eq 'EMS') ? 'EMS' : 'NAT';
}
$sth->finish();

# read domains positions from 
my %domains;
# if ($config->{domain_positions}) {
#   # get the positions of pfam domain annotations on transcripts
#   open (my $fh,"<",$config->{domain_positions});
#   while (<$fh>) {
#     chomp;
#     my ($tid,$start,$end,$domainType) = split /\t/, $_;
#     $domains{$tid} ||= [];
#     push @{$domains{$tid}}, {
#       start => $start,
#       end => $end,
#       type => $domainType
#     };
#   }
# }

$sth = $dbc->prepare(qq{
  SELECT vf.source_id,tv.feature_stable_id, tv.allele_string, tv.consequence_types, tv.cdna_start, sr.name, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, g.*
  FROM compressed_genotype_var g, variation_feature vf, transcript_variation tv, seq_region sr
  WHERE vf.seq_region_id = ?
  AND tv.consequence_types IN ($in_these_cs)
  AND vf.consequence_types NOT IN ($not_these_cs)
  AND tv.variation_feature_id = vf.variation_feature_id
  AND vf.seq_region_id = sr.seq_region_id
  AND vf.variation_id = g.variation_id
}, {'mysql_store_result' => 1});

$sth = $dbc->prepare(qq{
  SELECT vf.source_id,tv.feature_stable_id, tv.allele_string, tv.consequence_types, tv.cdna_start, sr.name, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, g.*
  FROM compressed_genotype_var g, variation_feature vf, transcript_variation tv, seq_region sr
  WHERE tv.consequence_types = ?
  AND tv.variation_feature_id = vf.variation_feature_id
  AND vf.seq_region_id = sr.seq_region_id
  AND vf.variation_id = g.variation_id
}, {'mysql_use_result' => 1});

my @json;
# foreach my $sr_id(@sr_ids) {
for my $con (@these_cs) {
    print STDERR "fetching $con variants\n";

#    $sth->execute($sr_id);
    $sth->execute($con);
  my $nrows = $sth->rows;
  # print STDERR "number of rows: $nrows\n";
  my %transcript_pop_consequence_sample;
  my ($source_id, $t, $a, $c, $cdna_start, $n, $s, $e, $r, $v, $ss, $g);
  $sth->bind_columns(\$source_id, \$t, \$a, \$c, \$cdna_start, \$n, \$s, \$e, \$r, \$v, \$ss, \$g);
  
  my %done;
  while($sth->fetch) {
      # print STDERR join ("\t", $source_id, $t, $a, $c, $cdna_start, $n, $s, $e, $r, $v, $ss),"\n";
    next unless $t2g{$t}; # only do canonical transcripts
    # my @consequences;
    # for my $con (split /,/, $c) {
    #     if (exists $consequence_level{$con} and $consequence_level{$con} > 0) {
    #         push @consequences, $con;
    #      }
    # }
    # print STDERR "2\n";
    # next unless @consequences;
    my ($ref,$alt,@etc) = split('/',$a);
    print STDERR "$a $ref, $alt, @etc\n" if (@etc);
	my @genotypes = unpack("(ww)*", $g);
	while(@genotypes) {
		my $sample_id = shift @genotypes;
		my $gt_code = shift @genotypes;
        # $sample_id or print STDERR "no sample id?\n";
        next unless $sample_id;
        next unless exists $genotype_to_ACGT{$gt_code}{$alt};
        my $ipop = $individual_pop{$sample_id};
        my $isource = $individual_source{$sample_id};
        # print STDERR "5\n";
        # next unless $ipop and $isource;
        # $isource or print STDERR "no isource for sample_id $sample_id\n";
        next unless $isource;
        # $isource == $source_id or print STDERR "source mismatch $sample_id, $isource, $source_id\n";
        # next unless $isource == $source_id;
        my $alleles = scalar keys %{$genotype_to_ACGT{$gt_code}} == 1 ? 'homo' : 'het';
        # print STDERR "7 $alleles @consequences\n";
        $transcript_pop_consequence_sample{$t2g{$t}}{$ipop}{$con}{$alleles}{$sample_id} = 1;
#        $transcript_pop_consequence_sample{$t2g{$t}}{$ipop}{$con}{$alleles}{$individual_id} = 1;
        # for my $con (@consequences) {
        #   if ($con eq 'missense_variant' and $domains{$t}) {
        #     my %dom_types;
        #     for my $dom (@{$domains{$t}}) {
        #       if ($dom->{start} <= $cdna_start and $cdna_start <= $dom->{end}) {
        #           $dom_types{$dom->{type}}=1;
        #         # $con .= "_".$dom->{type}
        #       }
        #     }
        #     if (keys %dom_types) {
        #         $con = join("_",$con, sort keys %dom_types);
        #     }
        #   }
        #   print STDERR "GOT $t2g{$t} $ipop $con $alleles $individual_id\n";
        #   $transcript_pop_consequence_sample{$t2g{$t}}{$ipop}{$con}{$alleles}{$individual_id} = 1;
        # }
    }
    # print join("\t", ($t2g{$t}, $t, $a, $c, $n, $s, $e, $r, join(',',@alt_samples))), "\n" if (@alt_samples);
  }
  $sth->finish();
  print STDERR "finished reading variants\n";
  for my $tid (sort keys %transcript_pop_consequence_sample) {
    my %any;
    my @veps;
    for my $pop_id (sort keys %{$transcript_pop_consequence_sample{$tid}}) {
      my $stype = $study_type{$pop_id};
      my $global_pop_id = $prod_name . "__" . $pop_id;
      for my $con (sort keys %{$transcript_pop_consequence_sample{$tid}{$pop_id}}) {
        for my $allele (keys %{$transcript_pop_consequence_sample{$tid}{$pop_id}{$con}}) {
          my @samples = sort keys %{$transcript_pop_consequence_sample{$tid}{$pop_id}{$con}{$allele}};
          # my @global_samples = map {$taxon_id . "." . $_} @samples;
          # push @veps, "\"vep__${global_pop_id}__${con}__${allele}\":[".join(',', @global_samples)."]\n";
          print join("\t", $tid, $global_pop_id, $stype, $con, $allele, join(",",map {$individual_name{$_}} @samples) ),"\n";
          # for my $sample (keys %{$transcript_pop_consequence_sample{$tid}{$pop_id}{$con}}) {
          #     print join("\t", $tid, $con, join("_",$taxon_id,$pop_id), $sample)),"\n";
            # $any{$pop_id}{$sample} = 1;
            # $any{any}{$sample} = 1;
         }
      }
    }
#    push @json, "\"$tid\":{\n".join(',',@veps)."}\n";
    # for my $pop_id (keys %any) {
    #   print ",\"vep__${pop_id}__any\":[",join(',', keys %{$any{$pop_id}}),"]";
    # }
    # print "}\n";
  }
} 
print STDERR "finished\n";
#print "{",join(",",@json),"}\n";
