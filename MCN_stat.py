# -*- coding: utf-8 -*-
from pymongo import MongoClient

client = MongoClient()
client = MongoClient('localhost', 27027)

dbName = "DatasetXmethods3"
db = client[dbName]


collection = "MCN_from_ETERNA_R95_0000_noPred-col-collectionName_stat_A_0.5_S_1.5_stn_1.5".replace(".","_").replace("-","_")
collection_out = collection +"_stat"



print ("connection")

array_so = db[collection].find({"soft":"so"}).distinct("ncm")

for ncm in array_so:
  soNcm_Low = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Low"})))
  soNcm_Bg = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Bg"})))
  soNcm_Hi = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Hi"})))
  print ("so : "+ncm)

  db[collection_out].insert({"ncm":ncm,"soft":"so","low":soNcm_Low,"bg":soNcm_Bg,"hi":soNcm_Hi})
  
  
  
array_mcff = db[collection].find({"soft":"mcff"}).distinct("ncm")

for ncm in array_mcff:
  mcffNcm_Low = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Low"})))
  mcffNcm_Bg = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Bg"})))
  mcffNcm_Hi = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Hi"})))
  print ("mcff : "+ncm)

  db[collection_out].insert({"ncm":ncm,"soft":"mcff","low":mcffNcm_Low,"bg":mcffNcm_Bg,"hi":mcffNcm_Hi})
  
print("fin : "+collection_out)