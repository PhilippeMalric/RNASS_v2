# -*- coding: utf-8 -*-

#ce script permet d'ecrire un fichier csv avec les information de la bas ede donner : ncm_best_stat

from pymongo import MongoClient


import sys
import os


from os import listdir
from os.path import isfile, join

import json

import statistics




'''
dirs = [o for o in os.listdir(d) ]
tab = []

for folderName in dirs:


  mypath = d+folderName

  #print("path : "+mypath) 
  if(os.path.isdir(mypath)):
    onlyfiles = [file for file in listdir(mypath+"/") if (isfile(join(mypath, file)) and file.endswith(".json"))]
    tab.append(len(onlyfiles))
    print("number of files : "+folderName+" : "+str(len(onlyfiles)))
  
print("sum : "+str(sum(tab)))
print("max : "+str(max(tab)))
print("min : "+str(min(tab)))
'''

client = MongoClient()
client = MongoClient('10.0.0.161', 27027)
dbName = "DatasetXmethods3"
db = client[dbName]
collection = "MCN_from_ETERNA_R95_0000_noPred-col-collectionName_stat_A_0.5_S_1.5_stn_1.5".replace(".","_").replace("-","_")




tab_so = []

array_so = db[collection].aggregate([{"$match":  {'soft': 'so'}},{"$group": {"_id": "null", "uniqueValues": {"$addToSet": "$ncm"}}}])


for ncm in array_so.next()["uniqueValues"]:
  soNcm = db[collection].find({"ncm":ncm,"soft":"so"})
  print ("so : "+str(soNcm[0]["low"]+soNcm[0]["bg"]+soNcm[0]["hi"]))

  tab_so.append([ncm,str((soNcm[0]["low"]+soNcm[0]["bg"]+soNcm[0]["hi"]))])
  
tab_mcff = []  
  
array_mcff = db[collection].aggregate([{"$match":  {'soft': 'mcff'}},{"$group": {"_id": "null", "uniqueValues": {"$addToSet": "$ncm"}}}])


for ncm in array_mcff.next()["uniqueValues"]:
  mcffNcm = db[collection].find({"ncm":ncm,"soft":"mcff"})
  print ("mcff : "+str(mcffNcm[0]["low"]+mcffNcm[0]["bg"]+mcffNcm[0]["hi"]))
  
  tab_mcff.append([ncm,str((mcffNcm[0]["low"]+mcffNcm[0]["bg"]+mcffNcm[0]["hi"]))])

d_out = '/u/malricp/bestRNA/csv/'
f = open(d_out+"ncmFreqDistribution_so.csv",'w')
f.write("ncm,freq\n")
f.write("\n".join(x[0]+","+x[1] for x in tab_so))
f.close()

f = open(d_out+"ncmFreqDistribution_mcff.csv",'w')
f.write("ncm,freq\n")
f.write("\n".join(x[0]+","+x[1] for x in tab_mcff))
f.close()


