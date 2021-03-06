// setup some dependencies
var express  = require('express'),
    compression = require('compression'),
    request = require('request'),
    qs = require('querystring'),
    cors = require('cors'),
    bodyParser = require('body-parser');

var app = express();

app.use(cors());
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended:true}));
app.use(compression());

var port = process.argv[2];
var solrURL = process.argv[3];

var cores = {
  genes : {
    description : "genes go here"
  },
  taxonomy : {
    description : "the ncbi taxonomy (partial)"
  },
  interpro : {
    description : "interpro domains"
  },
  GO : {
    description : "the gene ontology"
  },
  PO : {
    description : "the plant ontology"
  }
}

// define routes
app.get('/', function (req, res, next) { // return top level info
  res.json(cores);
});

app.get('/search/:core', function (req, res, next) {
  // check if the collection exists
  var c = req.params.core;
  if (!cores.hasOwnProperty(c)) res.json({"error":" '"+c+"' not found"});
  else {
    // consider validation of query params to secure solr...
    var url = solrURL
      + '/solr/' + c + '/query?'
      + qs.stringify(req.query);
    request.get(url).pipe(res);
  }
});

app.get('/suggest/:core', function (req, res, next) {
  // check if the collection exists
  var c = req.params.core;
  if (!cores.hasOwnProperty(c)) res.json({"error":" '"+c+"' not found"});
  else {
    // consider validation of query params to secure solr...
    var url = solrURL
      + '/solr/' + c + '/suggest?indent=true&wt=json&'
      + qs.stringify(req.query);
    request.get(url).pipe(res);
  }
});

app.post('/search/:core', function (req, res, next) {
  // check if the collection exists
  var c = req.params.core;
  if (!cores.hasOwnProperty(c)) res.json({"error":" '"+c+"' not found"});
  else {
    // consider validation of query params to secure solr...
    var url = solrURL
      + '/solr/' + c + '/query';
    req.body.wt = 'json';
    req.body.indent = true;
    request.post(url,{form:req.body}).pipe(res);
  }
});


var server = app.listen(port, function() {
  console.log('Listening on port %d', server.address().port);
});
