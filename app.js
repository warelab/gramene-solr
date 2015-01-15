// setup some dependencies
var express  = require('express'),
    request = require('request'),
    qs = require('querystring'),
    cors = require('cors'),
    bodyParser = require('body-parser'),
    validate = require('conform').validate;

// load settings from config file
var settings = require('./config/settings.json');

var app = express();

app.use(cors());
app.use(cors());
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended:true}));

var port = process.argv.length > 2 ? process.argv[2] : 8043;

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

app.get('/:core', function (req, res, next) {
  // check if the collection exists
  var c = req.params.core;
  if (!cores.hasOwnProperty(c)) res.json({"error":" '"+c+"' not found"});
  else {
    // consider validation of query params to secure solr...
    var url = settings.hostname
      + ':' + settings.port
      + '/solr/' + c + '/select?wt=json&indent=true&'
      + qs.stringify(req.query);
    request.get(url).pipe(res);
  }
});

app.post('/:core', function (req, res, next) {
  // check if the collection exists
  var c = req.params.core;
  if (!cores.hasOwnProperty(c)) res.json({"error":" '"+c+"' not found"});
  else {
    // consider validation of query params to secure solr...
    var url = settings.hostname
      + ':' + settings.port
      + '/solr/' + c + '/select';
    req.body.wt = 'json';
    req.body.indent = true;
    request.post(url,{form:req.body}).pipe(res);
  }
});


var server = app.listen(port, function() {
  console.log('Listening on port %d', server.address().port);
});
