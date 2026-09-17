> **Legacy manual recipe — superseded.** The build is driven by `build/` (`make solr`, stages
> `60_solr_genes` and `65_solr_suggestions`). The commands below are kept for reference only and are
> incomplete: they omit `TO.json`, `QTL_TO.json`, `qtls.json`, `TF.json`, `expr.json`, and the
> `Variants: rsID` layer (`suggestions/build_rsid_suggestions.sh`, ~10M docs, chunk-loaded because it
> is far too large for the single-curl loads shown here). They also predate the guards in the stages —
> core recreation, completeness gates, and the streaming loader. Use the stages.

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
