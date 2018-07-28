#!/usr/bin/env python
# -*- coding: utf-8 -*-
import statistics
import os
import RNASS
import json
import handler
from os import listdir
from os.path import isfile, join
import re

from pymongo import MongoClient
from optparse import OptionParser

from subprocess import Popen,PIPE

import pymongo

dir_of_prog = os.path.dirname(os.path.abspath(__file__)) +"/"

def scoreArrayAug(scoreA,index,f):
  if(len(scoreA) == 1):
    return [0]
  for i in range(0,int(index)-1):
    scoreA.insert(0,-999)
  while(len(scoreA) < f):
    scoreA.append(-999)  
  if(scoreA[-1] == 0):
    score = 0
    i = -1
    while(score == 0):
      scoreA[i] = -999
      i = i-1
      if(len(scoreA)>1 and i*-1 < len(scoreA)):
        score = scoreA[i]
      else: 
        return scoreA
  return scoreA[0:f]



#---------------------------for ED

 

def parseFold(lineSplited):
    finalTab = []
    
    freqMfe = -999
    ed = -999
    
    nextIsfreqMfe = False
    nextIsde = False
    for i in range(0,len(lineSplited)):
                #print("i : "+str(i)+" line : "+str(lineSplited[i]))
                if(nextIsfreqMfe):
                  freqMfe = float(lineSplited[i][:-1])
                  nextIsfreqMfe = False
                if(lineSplited[i] == 'ensemble' and freqMfe == -999):
                  nextIsfreqMfe = True
                
                if(nextIsde):
                  ed = float(lineSplited[i])
                  nextIsde = False
                if(lineSplited[i] == 'diversity'):
                  nextIsde = True
                
    return (ed,freqMfe)
       
       
       
def fastaWrite(outFile,seq,name):
  if not os.path.exists(outFile):
      try:  
              with open(outFile,"w") as newFasta:
                      header = ">" + name + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                      newFasta.write(header + seq)
      except IOError:
              print ("Failed to open " + outFile)
              exit(1)

def fastaErrase(outFile,name):
  if os.path.exists(outFile):
      try:  
              os.remove(outFile)
              
      except IOError:
              print ("Failed to errase : " + outFile)
              exit(1)

  psFile = dir_of_prog+"fasta_for_RNA_fold/"+name[:42]+"_ss.ps"
  if os.path.exists(psFile):
      try:  
              os.remove(psFile)
              
      except IOError:
              print ("Failed to errase : " + psFile)
              #exit(1)            
  else:
      print("remove failed : "+ psFile)

  #print ">> Done!"
def rnaFold(filePath,name,seq):
    os.chdir(filePath)
    fastaP = name+".fa"
    #print ("fastaP = "+fastaP)
    fastaWrite(fastaP,seq,name)
    stri = "RNAfold < \""+fastaP+"\" −noPS -p "
    #print (stri)
    Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
    output,err =  Process.communicate()
    fastaErrase(fastaP,name)
    os.chdir(dir_of_prog)
    #print ("out : "+output.strip().decode('ascii'))
    lineSplited = output.strip().decode('ascii').split()
    #la sequence
    #print lineSplited[3]
    #print ("line de Fold : "+str(len(lineSplited)))
    
    if(len(lineSplited)>0):
            seqFold = lineSplited[1]
            (ed,freqMfe) = parseFold(lineSplited[2:])
            
            #print("Ensemble diversity : "+str(ed)+" | frequence Mfe : "+str(freqMfe))
            
            #print("structureEnergieTab len : "+str(len(structureEnergieTab)))
    else:
      ed = 100000
      freqMfe = 100000
      print("no result for RNAfold")
    return (ed,freqMfe)

# arg parser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write report to FILE", metavar="FILE")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

ed_seuil = float(args[0])
stn_seuil = float(args[1])
score_seuil = float(args[2])

if(len(args) > 3):
    filesFilter = args[3]
else:
    filesFilter = "ETERNA_R00_0002"

if(len(args) > 4):
    collectionName = args[4]
else:
    collectionName = "OneEXPatAtime_Filter_ETERNA_R00_0002_noPred__A_0_5_S_5_0_stn_6_0_ed_5_0_stat"


prediction = True



keyWord = "S_pred"

dbName = "rdv"

if(prediction):
  p = re.compile('_S_(.*)_stn_(.*)_ed_(.*)')
  re_return = p.search(collectionName)
  pred_str = "_Pred_db-"+dbName+"-col-"+"_".join(re_return.group(x) for x in [1,2,3])+"_"
else:
  collectionName="-"
  pred_str = "_noPred_" 
  

adenineCutOff = 0.5
  
#filepath est le repertoire des fichier JSON  
id_unique = (keyWord+"_Filter_"+filesFilter+pred_str+"_A_"+str(adenineCutOff)+"_S_"+str(score_seuil)+"_stn_"+str(stn_seuil)+"_ed_"+str(ed_seuil)).replace(".","_").replace("-","_")
id_unique_short = keyWord+"_"+"_".join([str(score_seuil),str(stn_seuil),str(ed_seuil)])

if(options.verbose):
  print("options : "+str(options))
  print("args : "+str(args))

# fichier qui monitor certain parametre

f_rdat_associated_with_hi_stn = open(dir_of_prog+"meta_info/Hi_stn.txt","a")



'''
f = open("/u/malricp/_MCfoldVsRNAsubopt/projet/Mam4/folder2.txt",'r')
r = f.read()
f.close()
folders = r.split("\n")
'''

fa_adList_folder = dir_of_prog+"Dataset_1m7/"

'''
lineId = int(sys.argv[1])
lineId2 = int(sys.argv[2])
file = sys.argv[3]
'''
stn = 0
counter =0

rna2DTab = []
folder = "/u/malricp/_MCfoldVsRNAsubopt/projet/Mam4/rmdb/"
onlyfiles = [f for f in listdir(folder) if isfile(join(folder, f))]

folder_in = "rmdb"
linesParsed = []
for file in onlyfiles :
  if(options.verbose):
    print("fileFilter : "+filesFilter + " file : "+file)
  if(file.startswith(filesFilter)):
      root = file[:-5]
      if(options.verbose):
        print("root : "+root )
      rdat = handler.RDATFile()
      rdat.load(open(folder+file))
      rdat.validate()
      global_annotation = rdat._annotation_str(rdat.annotations,"|")
      for (key,value) in rdat.constructs.items() :
        if(hasattr(rdat.constructs[key].data[0],"annotations")):
          if(options.verbose):
            print("annotation in data")
            print("annotations : "+str(rdat.constructs[key].data[0].annotations))
          annot_tab_str = [x.annotations for x in rdat.constructs[key].data]
          if("sequence" in rdat.constructs[key].data[0].annotations):
              if(options.verbose):
                print ("sequence in annotation")
              rdat.seqTab = [x.annotations["sequence"][0] for x in rdat.constructs[key].data]
              #print("seqTab : "+"\n".join(rdat.seqTab))
              if(options.verbose):
                print("seqTab len : "+str(len(rdat.seqTab)))
              for i in range(0,len(rdat.seqTab)):
                seq = rdat.seqTab[i]
                (ed,freqMfe) = rnaFold(dir_of_prog+"fasta_for_RNA_fold/",root+"_"+str(i)+id_unique,seq)
                v = value.data[i].values
                v_avg = statistics.mean(v)
                e = value.data[i].errors
                offS = rdat.constructs[key].offset
                v = scoreArrayAug(v,offS,len(seq))
                e = scoreArrayAug(e,offS,len(seq))
                if(options.verbose):
                  print("ed : "+str(ed)+" v_avg : "+str(v_avg))
                if("signal_to_noise" in rdat.constructs[key].data[0].annotations):
                  if(options.verbose):
                    print ("signal_to_noise in annotation")
                  stn = rdat.constructs[key].data[i].annotations["signal_to_noise"][0].split(":")[1]
                else : 
                  if(options.verbose):
                    print ("No signal_to_noise in annotation")
                  stn = 100 
                  
                c1 = ("modifier" in rdat.constructs[key].data[0].annotations and rdat.constructs[key].data[i].annotations["modifier"][0] == "1M7")
                c2 = ("modifier" in rdat.annotations and rdat.annotations["modifier"][0] == "1M7")
                
                if(float(stn) < stn_seuil and v_avg < score_seuil and ed < ed_seuil):
                  if(c1 or c2 ): 
                    if(options.verbose):
                      print("1m7")
                    oneLineParsed =  {}
                    oneLineParsed.clear()
                    oneLineParsed["ed"] = ed   
                    oneLineParsed["stn"] = stn   
                    oneLineParsed["freqMfe"] = freqMfe   
                    oneLineParsed["seq"] = seq
                    oneLineParsed["structure"] = ""
                    oneLineParsed["scoreTab"] =  v
                    oneLineParsed["errorTab"] = e
                    oneLineParsed["old_TSV_ID"] = "id_"+str(counter)
                    oneLineParsed["old_info"] = rdat.constructs[key].data[i].annotations
                    #print("stn : "+str(rnaOld["stn"]))
                    linesParsed.append(oneLineParsed)
                    counter += 1
                    if(counter > 1000000):
                      print("break")
                      break
                else:
                  f_rdat_associated_with_hi_stn.write("file2_"+root+"__stn:"+stn+"\n")
                  if(options.verbose):
                    print("file : "+root+"__stn: "+str(stn)+"\n")
      #-----------------------------------------------------------------------------
      if(len(linesParsed)>0):
        
        nbNt_Ncm = 10


        filePath = dir_of_prog+"Dataset_1m7/"+id_unique+"/"
        if(options.verbose):
          print("filePath : "+filePath)
        if (not os.path.exists(filePath)):
          os.makedirs(filePath)
        #RNA_featureChoice.RNA2D(filePath,root,linesParsed,folder,"experience_6_nov_2017",id_num) 
        d = {
            'inputFileLoader': "test",
            'filePath': fa_adList_folder,
            'publicFolder': filePath,
            'currentExp': "experience_24_juillet_2018",
            'uniqueId': "0",
            'soTreshold': 1,
            'mcffTreshold':1,
            'CSVFolder':"CSV_FOLDER/"+"test",
            'prediction':prediction ,
            'printCsv':False,
            'so_e_value':1,
            'lowScore' : 0.5,
            'hiScore' : 1,
            'predCutOff' : 0.5,
            'stnSeuil' : stn_seuil,
            'nbNt_Ncm' : nbNt_Ncm,
            "adenineCutOff" : adenineCutOff,
            "scoreCutOff" : score_seuil,
            "ed_seuil" : ed_seuil,
            "collectionName" : collectionName,
            "dbName":dbName,
            "root":root,
            "MFP":False
            }
        if(options.verbose):
          print("linesParsed : "+str(len(linesParsed)))
        RNASS.RNASS(d,linesParsed,options)
        if(options.verbose):
          print(id_unique)
      else:
        print("noSeq")


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

if (prediction):
    if(options.verbose):
      print("startLoaderCompile")
    client = MongoClient()
    client = MongoClient('majsrv1', 27027)
    dbName = "rdv"
    db = client[dbName]
    if(options.verbose):
      print("connected")

    collection = id_unique


    publicFolder = dir_of_prog+"Dataset_1m7/"+id_unique

    pathTab = []
    pathTab.append({"name":"so_detail_mfe","path":["ncmTabDG_so"],"soft":"so"})
    pathTab.append({"name":"mcff_detail_mfe","path":["ncmTabDG_mcff"],"soft":"mcff"})
    #pathTab.append({"name":"so_ncm","path":["localNcmD","so"],"soft":"so"})
    #pathTab.append({"name":"mcff_ncm","path":["localNcmD","mcff"],"soft":"mcff"})


    onlyfiles = [f for f in listdir(publicFolder)]

    #onlyfiles = [file for file in listdir(mypath) if isfile(join(mypath, file))]

    if(options.verbose):
      print("onlyfiles : "+ str(len(onlyfiles)))

    rnaTab = []
    counter = 0

    for file in onlyfiles:
        if(file.endswith(".json") ):
          counter += 1
          print("file i : "+str(counter))
          if(counter > 10000000):
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


    if(options.verbose):
      print ("connection")

    array_so = db[collection].find({"soft":"so"}).distinct("ncm")
    counter = 0
    for ncm in array_so:
      counter += 1
      print(str(counter)+") ncm : "+ncm)
      soNcm_Low = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Low"})))
      soNcm_Bg = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Bg"})))
      soNcm_Hi = len(list(db[collection].find({"ncm":ncm,"soft":"so","label":"Hi"})))
      #print ("so : "+ncm)

      db[collection_out].insert({"ncm":ncm,"soft":"so","low":soNcm_Low,"bg":soNcm_Bg,"hi":soNcm_Hi})



    array_mcff = db[collection].find({"soft":"mcff"}).distinct("ncm")
    counter = 0
    for ncm in array_mcff:
      counter += 1
      print(str(counter)+") mcff_ncm : "+ncm)
      mcffNcm_Low = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Low"})))
      mcffNcm_Bg = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Bg"})))
      mcffNcm_Hi = len(list(db[collection].find({"ncm":ncm,"soft":"mcff","label":"Hi"})))
      #print ("mcff : "+ncm)

      db[collection_out].insert({"ncm":ncm,"soft":"mcff","low":mcffNcm_Low,"bg":mcffNcm_Bg,"hi":mcffNcm_Hi})

    db[collection_out].create_index(
        [("soft",pymongo.DESCENDING),("ncm",pymongo.DESCENDING)],
        unique=True
    )

f_rdat_associated_with_hi_stn.close()
print("fin")     
print("collection : " + collection_out)




print("options : "+str(options))

print("options verbose : "+str(options.verbose))

print("args : "+str(args))
