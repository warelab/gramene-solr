
// Your JSON data
const data = require('./qtls.json');

// Solr core URL
const solrCoreUrl = 'http://squam:8983/solr/sorghum_genes9';

// Function to create a query string
const createQueryString = (params) => {
  return Object.entries(params)
    .map(([key, value]) => `${encodeURIComponent(key)}=${encodeURIComponent(value)}`)
    .join('&');
};

// Function to query Solr
const querySolr = async (fqField, fqValue) => {
  try {
    const params = {
      q: '*:*', // General query to match all documents
      fq: `${fqField}:${fqValue}`,
      fl: 'taxon_id',
      rows: 1,  // Limit to 1 document
      wt: 'json' // JSON response
    };
    const queryString = createQueryString(params);
    
    const response = await fetch(`${solrCoreUrl}/select?${queryString}`);
    if (!response.ok) {
      throw new Error(`HTTP error! Status: ${response.status}`);
    }

    const responseData = await response.json();
    const numFound = responseData.response.numFound;
    const firstDoc = responseData.response.docs.length > 0 ? responseData.response.docs[0] : null;

    return { numFound, firstDoc };
  } catch (error) {
    console.error('Error querying Solr:', error.message);
    return { numFound: 0, firstDoc: null };
  }
};

// Iterate over the JSON objects and query Solr
const processQTLs = async () => {
  let qtls = [];
  for (const obj of data) {
    console.error(`Processing ID: ${obj.id}`);
    const result = await querySolr(obj.fq_field, obj.fq_value);

    console.error(`ID: ${obj.id}, Total Results: ${result.numFound}`);
    obj.num_genes = result.numFound;
    if (result.firstDoc) {
      console.error(`First Document: ${JSON.stringify(result.firstDoc, null, 2)}`);
      obj.taxon_id = [result.firstDoc.taxon_id];
      obj.taxon_freq = [result.numFound];
      qtls.push(obj);
      // console.log(JSON.stringify(obj));
    } else {
      console.error('No documents found for this query.');
    }
  }
  console.log(JSON.stringify(qtls,null,2));
};

// Run the script
processQTLs();
