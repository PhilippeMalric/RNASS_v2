//-------------------------------------------------api
var express = require("express");
var app = express()
var path = require("path");
var bodyParser = require("body-parser");
var mongodb = require("mongodb");
var cors = require('cors')

RNAS_COLLECTION = "RNA"

app.use(cors())

app.use(bodyParser.json());

var db;


// Connect to the database before starting the application server.
mongodb.MongoClient.connect("mongodb://localhost:27027/rdv", function (err, client) {
  if (err) {
    console.log(err);
    process.exit(1);
  }

  // Save database object from the callback for reuse.
  db = client.db('rdv');
  console.log("Database connection ready");

  // Initialize the app.
  var server = app.listen(3000, function () {
    var port = server.address().port;
    console.log("App now running on port", port);
  });
});

app.get("/test", function (req, res) {
  res.send("test")
})


app.get("/ncm_grouped_Low_std_dev_so/collection=:collections/skip=:skip/limit=:limit/countMin=:cmin/stdDevMax=:stdDevMax", function (req, res) {
  //console.log("id :" + req.params.id)

  db.collection(req.params.collections).aggregate([{ $project: { nts: 1, info:{root:1} } },
  { $unwind: { path: "$nts" } },
  { $project: { nts: { $arrayElemAt: ["$nts.ncmTabDG_so", 0] }, score: "$nts.score", rna: "$nts.rna_id",pos: "$nts.position",exp:"$info.root" } },
  { $project: { ncm: "$nts.ncm.merge", score: 1,rna:1,pos:1,exp:1} },
  {
    $group: {
      _id: "$ncm",
      stdDev: { $stdDevPop: "$score" },
      scoreMoy: { $avg: "$score" },
      count: { $sum: 1 },
      items: { $push:  { score: "$score", rna: "$rna",pos:"$pos",exp:"$exp" }}
    }
  },
  {
    '$match': {
      $and: [
        { _id: { $nin: [/.*L.*/i] } },
        { 'stdDev': { '$lt': Number(req.params.stdDevMax) } },
        { 'stdDev': { '$ne': 0 } },
        { 'count': { '$gt': Number(req.params.cmin) } }]
    }
  },
  { $skip: Number(req.params.skip) },
  { $limit: Number(req.params.limit) }])
    .toArray(function (err, docs) {
      if (err) {
        handleError(res, err.message, "Failed to get contacts.");
      } else {
        res.status(200).json(docs);
      }
    });
});

app.get("/ncm_grouped_Low_std_dev_mcff/collection=:collections/skip=:skip/limit=:limit/countMin=:cmin/stdDevMax=:stdDevMax", function (req, res) {
  //console.log("id :" + req.params.id)

  db.collection(req.params.collections).aggregate([{ $project: { nts: 1, info:{root:1} } },
  { $unwind: { path: "$nts" } },
  { $project: { nts: { $arrayElemAt: ["$nts.ncmTabDG_mcff", 0] }, score: "$nts.score", rna: "$nts.rna_id",pos: "$nts.position",exp:"$info.root" } },
  { $project: { ncm: "$nts.ncm.merge", score: 1,rna:1,pos:1,exp:1} },
  {
    $group: {
      _id: "$ncm",
      stdDev: { $stdDevPop: "$score" },
      scoreMoy: { $avg: "$score" },
      count: { $sum: 1 },
      items: { $push:  { score: "$score", rna: "$rna",pos:"$pos",exp:"$exp" }}
    }
  },
  {
    '$match': {
      $and: [
        { _id: { $nin: [/.*L.*/i] } },
        { 'stdDev': { '$lt': Number(req.params.stdDevMax) } },
        { 'stdDev': { '$ne': 0 } },
        { 'count': { '$gt': Number(req.params.cmin) } }]
    }
  },
  { $skip: Number(req.params.skip) },
  { $limit: Number(req.params.limit) }])
    .toArray(function (err, docs) {
      if (err) {
        handleError(res, err.message, "Failed to get contacts.");
      } else {
        res.status(200).json(docs);
      }
    });
});






app.get("/ncm_stat/collection=:collection/soft=:soft/minimum=:min", function (req, res) {
  //console.log("id :" + req.params.id)

  db.collection(req.params.collection).aggregate([{ $project: { "ncm":1,"soft":1,"low":1,"bg":1,"hi":1 } },
  {
   '$match':
    {$and: [
     {soft:req.params.soft},
     {$or :[{ 'hi': { '$gt': Number(req.params.min) } },
            { 'low': { '$gt': Number(req.params.min) } }]}
      ]}}])
    .toArray(function (err, docs) {
      if (err) {
        handleError(res, err.message, "Failed to get contacts.");
      } else {
        res.status(200).json(docs);
      }
    });
});


app.get("/rna/collection=:collection/exp=:exp/id=:id", function (req, res) {
  //console.log("id :" + req.params.id)

  db.collection(req.params.collection).find({rna_id:req.params.id,"info.root":eq.params.exp})
    .toArray(function (err, docs) {
      if (err) {
        handleError(res, err.message, "Failed to get contacts.");
      } else {
        res.status(200).json(docs);
      }
    });
});










