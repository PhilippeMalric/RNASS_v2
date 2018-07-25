# -*- coding: utf-8 -*-


import sys
import os


from os import listdir
from os.path import isfile, join

import json

import statistics

from pymongo import MongoClient

def filterMinus999(tab):
  newTab = []
  for x in tab:
    if(x != -999):
      newTab.append(x)
  return newTab

def filterDash(tab):
  tabtemp = []
  for x in tab:
    if(x[0] != "-"):
      tabtemp.append(x)
  return tabtemp


print("startLoaderCompile")
client = MongoClient()
client = MongoClient('localhost', 27027)
dbName = "rdv"
db = client[dbName]
collection = "ETERNA_R95_0000_Pred_db_dbName_col_collectionName_stat_A_0_5_S_1_5_stn_5".replace(".","_").replace("-","_")
print("connected")

publicFolder = "/u/malricp/DatasetXmethods/Dataset_1m7/ETERNA_R95_0000_Pred_db-dbName-col-collectionName_stat_A_0.5_S_1.5_stn_5"

pathTab = []
pathTab.append({"name":"so_detail_mfe","path":["ncmTabDG_so"],"soft":"so"})
pathTab.append({"name":"mcff_detail_mfe","path":["ncmTabDG_mcff"],"soft":"mcff"})
#pathTab.append({"name":"so_ncm","path":["localNcmD","so"],"soft":"so"})
#pathTab.append({"name":"mcff_ncm","path":["localNcmD","mcff"],"soft":"mcff"})


onlyfiles = [f for f in listdir(publicFolder)]

#onlyfiles = [file for file in listdir(mypath) if isfile(join(mypath, file))]

print("onlyfiles : "+ str(len(onlyfiles)))

rnaTab = []
counter = 0

for file in onlyfiles:
    if(file.endswith(".json") ):
      print("i : "+str(counter))
      counter += 1
      if(counter > 100000):
        break
      rna = json.loads(open(publicFolder+"/"+file,"r").read())
      for dPath in pathTab:
        d = {}
        rnaD = {}
        t = rna["nts"]
        filtered = filterMinus999(rna["scoreTab"])
        hi_threshold = 1
        low_threshold = 0.5 
        for nt in t:
          o = nt
          for p in dPath["path"]:
            #print("keys : "+str(o.keys()))
            mcn = o[p]
          url = "http://majsrv1.iric.ca:3000/"+"|"+str(rna["rna_id"])+"|"+str(nt["position"])
          #tbTemp.append(str(v))
          mcn = mcn[0]["ncm"]["merge"]
          root = mcn.replace(" ","&")
          freq = "mfe"
          score = str(nt["score"])
          so_pairee = str(nt["freqPairee_so"])
          mcff_pairee = str(nt["freqPairee_mcff"])
          if(nt["score"] < low_threshold and nt["score"] != -999 ):
            label ="Low"
            db[collection].insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
          
          if(nt["score"] > hi_threshold and nt["score"] != -999):
            label = "Hi"
            db[collection].insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
            
          if((nt["score"] < hi_threshold and nt["score"] != -999) and nt["score"] > low_threshold ):
            label = "Bg"
            db[collection].insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
          
          
          
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
  
print("fin")          
   
  