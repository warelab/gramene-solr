convert mongo genes to solr genes
```
cd genes
node --max-old-space-size=8192 ./mongo2solr.js > solr_genes.json
```

load into solr genes core
```
curl 'http://brie:8984/solr/genes55/update?commit=true' --data-binary @solr_genes.json -H 'Content-type:application/json'
```

create suggestions
```
cd ../suggestions
node aux.js http://brie:8984/solr/genes55
node --max-old-space-size=8192 ./genes.js http://brie:8984/solr/genes55 > genes.json
```

load into solr suggestions core
```
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @GO.json -H 'Content-type:application/json'
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @PO.json -H 'Content-type:application/json'
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @taxonomy.json -H 'Content-type:application/json'
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @domains.json -H 'Content-type:application/json'
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @pathways.json -H 'Content-type:application/json'
curl 'http://brie:8984/solr/suggestions55/update?commit=true' --data-binary @genes.json -H 'Content-type:application/json'
```
