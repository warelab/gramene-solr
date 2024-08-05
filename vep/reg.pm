use Bio::EnsEMBL::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);


my $def_user = 'weix';
my $def_pass = 'warelab';
my $def_host = 'cabot';
my $def_port = 3306;


Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'sorghum_bicolor',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'sorghum_bicolor_core_8_108_30'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'sorghum_bicolor',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'sorghum_bicolor_variation_8_108_30'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa_core_8_108_7'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa_variation_8_108_7'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_aus',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_aus_core_8_108_2'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_aus',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_aus_variation_8_108_2'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa117425',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa117425_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa117425',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa117425_variation_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa125827',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa125827_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativa125827',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativa125827_variation_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativaazucena',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativaazucena_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativaazucena',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativaazucena_variation_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativair64',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativair64_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativair64',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativair64_variation_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativamh63',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativamh63_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativamh63',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativamh63_variation_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativazs97',
	-group   => 'core',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativazs97_core_8_108_1'
);

Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-species => 'oryza_sativazs97',
	-group   => 'variation',
	-port    => $def_port,
	-host    => $def_host,
	-user    => $def_user,
	-pass    => $def_pass,
	-dbname  => 'oryza_sativazs97_variation_8_108_1'
);


1;



