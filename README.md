first dump the genes from mongo (or just use the genes that you already have locally)
```
mongoexport -c genes > mongo_genes.json
```

convert mongo genes to solr genes
```
cd genes
./mongo2solr.js > solr_genes.json
```

load into solr genes core
```
curl 'http://brie:8983/solr/genes50/update?commit=true' --data-binary @solr_genes.json -H 'Content-type:application/json'
```

create suggestions
```
cd ../suggestions
node aux.js http://brie:3983/solr/genes50
mongoexport -h brie -d search50 -c genes | --max-old-space-size=8192 node genes.js http://brie:8983/solr/genes50 > genes.json
```

load into solr suggestions core
```
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @GO.json -H 'Content-type:application/json'
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @PO.json -H 'Content-type:application/json'
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @taxonomy.json -H 'Content-type:application/json'
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @domains.json -H 'Content-type:application/json'
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @pathways.json -H 'Content-type:application/json'
curl 'http://brie:8983/solr/suggestions50/update?commit=true' --data-binary @genes.json -H 'Content-type:application/json'
```
