first dump the genes from mongo (or just use the genes that you already have locally)
```
mongoexport -c genes > mongo_genes.json
```

convert mongo genes to solr genes
```
cd genes
cat mongo_genes.json | ./mongo2solr.js > solr_genes.json
```

load into solr genes core
```
curl 'http://localhost:8983/solr/genes/update?commit=true' --data-binary @solr_genes.json -H 'Content-type:application/json'
```

create suggestions
```
cd ../suggestions
mongoexport -h brie -d search48 -c GO | ./aux.js GO http://localhost:8983/solr/genes' > GO.json
mongoexport -h brie -d search48 -c PO | ./aux.js PO http://localhost:8983/solr/genes' > PO.json
mongoexport -h brie -d search48 -c taxonomy | ./aux.js taxonomy http://localhost:8983/solr/genes' > taxonomy.json
mongoexport -h brie -d search48 -c domains | ./aux.js domains http://localhost:8983/solr/genes' > domains.json
cat mongo_genes.json | node genes.js http://localhost:8983/solr/genes > genes.json
```

load into solr suggestions core
```
curl 'http://localhost:8983/solr/suggestions48/update?commit=true' --data-binary @GO.json -H 'Content-type:application/json'
curl 'http://localhost:8983/solr/suggestions48/update?commit=true' --data-binary @PO.json -H 'Content-type:application/json'
curl 'http://localhost:8983/solr/suggestions48/update?commit=true' --data-binary @taxonomy.json -H 'Content-type:application/json'
curl 'http://localhost:8983/solr/suggestions48/update?commit=true' --data-binary @domains.json -H 'Content-type:application/json'
curl 'http://localhost:8983/solr/suggestions48/update?commit=true' --data-binary @genes.json -H 'Content-type:application/json'
```
