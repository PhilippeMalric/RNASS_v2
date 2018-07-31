# -*- coding: utf-8 -*-
import statistics
import os
import itertools
import re
import math
import json
import itertools
from subprocess import Popen,PIPE
from pymongo import MongoClient
import timeit

start = timeit.default_timer()


#from scipy import stats

#from Graph import Graph

#from Score import Score


#ecrit des fichier csv avec des choix de caracteristique different
# -2 - position relative
# -1 - sorte seulement
# 0 - subseq only
# 0.1 - pairee ou non
# 1 - ncm only
# 2 - in helice only
# 3 - reduce set of ncmCaracteristic_d
# 4 - All in


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class RNASS:
      'generate an html file with d3.js in it to visualize reactivity data on secondary structure'
      count = 0 
      

#filePath,name,seqProfilTab,folder,publicF,expName,id_num
      
      def __init__(self,d,linesParsed,options):
          self.options = options
          if(options.verbose):
            print("options.verbose : "+str(options.verbose))
          
          self.fileError = open("fileError.txt",'a')
          
          if(options.verbose):
            print ("begining of RNA.py -----------------------------")
          RNASS.count += 1  
          if(options.verbose):
            print("RNASS.count : "+str(RNASS.count ))
            print(str(d))
          #populate self variable----------------------------------------------
          self.defVar(d,linesParsed)
          
    
            
          # pour determiner si un nucleotide est la majorite du temps dans une helice
          InH_sens = 0.75

          #unique pour les ARN et pour les nucleotides
          rna_id_num = -1
          nt_id_num = -1
          
          if(options.verbose):
            print("self.seqProfilTab len : " + str(len(self.seqProfilTab)))
          
          for e in self.seqProfilTab:
                  rna_id_num += 1
                  if(options.verbose):
                    print(e["seq"])
                    print("############################ rna_id_num : "+str(rna_id_num))
                  stop1 = timeit.default_timer()
                  if(options.verbose):
                    print (str(stop1 - start))
                  #Your statements here


                  self.seq = e["seq"]
                  #print("seqUnderTreatment : "+e["seq"])
                  
                  #Un id unique pour chaque ARN d'un fichier rdat
                  self.id_ARN = int(self.id_ARN)+ 1
                  #print("id_ARN : "+str(self.id_ARN))
                  
                  #pour enlever les nucleotides qui n'ont pas de score
                  filtered = self.filterMinus999(e["scoreTab"])
                  if(len(filtered)>1): 
                    self.nts_score_mean =  statistics.mean(filtered)
                    self.nts_score_sd = statistics.stdev(filtered)
                  else:
                    if(options.verbose):
                      print("scoreTab : "+str(e["scoreTab"]))
                    continue
                  #filtre 
                  if((float(e["stn"]) > float(self.d['stnSeuil'])) or (float(self.nts_score_mean) > float(self.scoreCutOff)) or (float(self.calculFreqR("A",e["seq"])) > self.adenineCutOff)):
                    if(options.verbose):
                      print("freqA : "+str(self.calculFreqR("A",e["seq"])))
                      print("self.nts_score_mean : "+str(self.nts_score_mean))
                      print("stn : "+e["stn"])
                    continue
                  if(options.verbose):
                    print("############################ Apres le filtre sur les Adenine ("+str(self.calculFreqR("A",e["seq"]))+") et le score "+str(self.nts_score_mean) +"et le stn : "+str(e["stn"]))
                  stop2 = timeit.default_timer()
                  if(options.verbose):
                    print (str(stop2 - start))
                  #print("e['scoreTab'] : "+str(e["scoreTab"]))   
                  self.reactivityVector = []
                  for sc in e["scoreTab"]:
                      #print("sc : "+str(sc))
                      #print("lowScore : "+str(lowScore))
                      if(sc > self.d['hiScore']):
                          self.reactivityVector.append("Hi")
                      else: 
                          if(sc < self.d["lowScore"] and sc != -999):
                              self.reactivityVector.append("Low")
                          else:
                              if(sc != -999):
                                  self.reactivityVector.append("Bg")
                              else:
                                  self.reactivityVector.append("NA")
                                  
                  if(options.verbose):
                    print("########## RNASO prediction")
                  # lance RNAsubopt et parse le resultat, une liste de structure/energie de la classe strucEnergy est retournezzz
                  self.sESTabsubOpt = []
                  self.sESTabsubOpt = self.rnaSO(e["seq"],True)
                  if(options.verbose):
                    print("self.sESTabsubOpt : "+ str(len(self.sESTabsubOpt)))
                  self.sESTabMcff = []
                  self.sESTabMcff = self.mcff(e["seq"],True) 
                  if(options.verbose):
                    print("self.sESTabMcff : "+ str(len(self.sESTabMcff)))
                  
                  
                  #Donne la frequence de pairage des nt
                  frP_mcff = self.sumNtWise([x.struct for x in self.sESTabMcff]) 
                  frP_so = self.sumNtWise([x.struct for x in self.sESTabsubOpt])
                  
                  # si le parsing des deux logiciel de prediction de la SS est complet et reussi : 
                  if(len(self.sESTabsubOpt) > 0 and len(self.sESTabMcff) > 0):
                      self.pairTabAllSubOpt_so = []
                      self.pairTabAllSubOpt_so = [self.findPairGen(x.struct) for x in self.sESTabsubOpt]
                      jsNode_so = [self.jsNodeAll(pairTab,len(e["scoreTab"])) for pairTab in self.pairTabAllSubOpt_so]
                      #print("jsNodePrint_so: "+jsNodePrint_so)
                      self.pairTabAllSubOpt_mcff = []
                      self.pairTabAllSubOpt_mcff = [self.findPairGen(x.struct) for x in self.sESTabMcff]
                      jsNode_mcff = [self.jsNodeAll(pairTab,len(e["scoreTab"])) for pairTab in self.pairTabAllSubOpt_mcff]
                      #print("jsNodePrint_mcff: "+jsNodePrint_mcff)
                      
                      if(options.verbose):
                        print("############################ Apres le Avoir trouve les pair") 
                      stop3 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop3 - start))
                      
                      if(self.MFP):
                        (t1_mcff,mfp_mcff) = self.findMostfreqPartner(self.pairTabAllSubOpt_mcff,len(self.sESTabMcff) )
                        (t1_so,mfp_so) = self.findMostfreqPartner(self.pairTabAllSubOpt_so,len(self.sESTabsubOpt))
                      else:
                        mfp_mcff = [0]*len(self.seq)
                        mfp_so = [0]*len(self.seq)
                      if(options.verbose):
                        print("############################ Apres le Avoir trouve les partenaires frequents") 
                      stop4 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop4 - start))
                      #Les predictionS se font ici
                      #--------------------------------------------------------------
                      
                      ncm_by_nt_Dict = {}
                      ncm_by_nt_Dict.clear()
                      ncm_by_nt_Dict = self.pbToNcmAll(False,False)
                      
                      ncm_by_nt_mfe_Dict = {}
                      ncm_by_nt_mfe_Dict.clear()
                      ncm_by_nt_mfe_Dict = self.pbToNcmAll(False,True)
                      
                      ncm_by_nt_d_Dict = {}
                      ncm_by_nt_d_Dict.clear()
                      ncm_by_nt_d_Dict = self.pbToNcmAll(True,False)
                      
                      ncm_by_nt_d_mfe_Dict = {}
                      ncm_by_nt_d_mfe_Dict.clear()
                      ncm_by_nt_d_mfe_Dict = self.pbToNcmAll(True,True)
                      #print("nt[isENANtsV] : "+ ",".join([str(x) for x in ncm_by_nt_d_Dict["isENANtsV"]]))
                      
                      if(options.verbose):
                        print("############################ Apres le Avoir trouve les ncm") 
                      stop5 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop5 - start))
                      #-----------------------------------------------------------
                      
                      
                      #print("tabVoisinAllsub_so[10] : "+str(self.tabVoisinAllsub_so[10]))
                      #print("tabVoisinAllsub_mcff[10] : "+str(self.tabVoisinAllsub_mcff[10]))
                      
                      dotBraquetMcffTab = [x.struct for x in self.sESTabMcff]
                      dotBraquetSoTab = [x.struct for x in self.sESTabsubOpt ]
                      
                      seq = e["seq"]
                      
                      
                      dotBr2p = self.jsDictGeneratorByPredictor(dotBraquetMcffTab,dotBraquetSoTab)
                      d3ForceLayout2p = self.jsDictGeneratorByPredictor(jsNode_mcff,jsNode_so)
                      
                      
                      _2NtMotif = self.extractMotif(2,seq)
                      #_2NtMoti#file = "Mapseek1.txt"f_str = json.dumps(_2NtMotif)
                      scoreTabStr = "["+",".join(['%.5f' % x for x in e["scoreTab"]])+"]"
                      errorTabStr = "["+",".join(['%.5f' % x  for x in e["errorTab"]])+"]"

                      # pour tous les nucleotide -------------------------------
                      nTs = self.nT_Treatment(seq,ncm_by_nt_d_Dict,e,frP_mcff,frP_so,nt_id_num,mfp_mcff,mfp_so)
                      
                      if(options.verbose):
                        print("############################ Apres le traitement des nucleotides") 
                      stop6 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop6 - start))
                      scoreT = [x["score"] for x in nTs]
                      
                      filtered = self.filterMinus999(scoreT)
                      if(len(filtered)>1): 
                        self.nts_score_mean =  statistics.mean(filtered)
                        self.nts_score_sd = statistics.stdev(filtered)
                      else:
                        if(options.verbose):
                          print("score T  : "+str(scoreT))
                        break
                      
                          
                      
                       
                       
                      # CONVERTI LA MATRICE D'ADJACENCE en liste de noeud et de lien pour le layout force de d3js
                      graph_transition = {}
                      (graph_transition["nodes"],
                       graph_transition["links"],
                       graph_transition["maxMcff"],
                       graph_transition["minMcff"],
                       graph_transition["maxSO"],
                       graph_transition["minSO"]) = self.mAdToD3Graph(
                                                                  self.sESTabMcff,
                                                                  self.sESTabsubOpt,
                                                                  self.adjListMcff,
                                                                  self.adjListSO)
                      
                      
                      self.RNA["oldInfo"] = e["old_info"]
                      self.RNA["info"] = self.d
                      self.RNA["rna_id"] = self.id_ARN
                      self.RNA["nts"] = nTs
                      self.RNA["ScoresMoy"] = self.nts_score_mean
                      self.RNA["frequenceNT"] = {}
                      self.RNA["frequenceNT"]["A"] = self.calculFreqR("A",seq)
                      self.RNA["frequenceNT"]["C"] = self.calculFreqR("C",seq)
                      self.RNA["frequenceNT"]["G"] = self.calculFreqR("G",seq)
                      self.RNA["frequenceNT"]["U"] = self.calculFreqR("U",seq)
                      if(self.MFP):
                        self.RNA["t1_mcff"] = t1_mcff
                        self.RNA["t1_so"] = t1_so
                      self.RNA["structure"] = e["structure"]
                      self.RNA["stn"] = e["stn"]
                      self.RNA["frP_mcff"] = frP_mcff
                      self.RNA["frP_so"] = frP_so
                      self.RNA["scoreTab"] = e["scoreTab"]
                      self.RNA["erreurTab"] = e["errorTab"]
                      self.RNA["seq"] = e["seq"]
                      self.RNA["ed"] = e["ed"]
                      self.RNA["freqMfe"] = e["freqMfe"]
                      self.RNA["graph_transition"] = graph_transition
                      self.RNA["dotBr2p"] = dotBr2p
                      self.RNA["d3ForceLayout2p"] = d3ForceLayout2p
                      self.RNA["_2NtMotif"] = _2NtMotif
                      self.RNA["scoreEcart_type"] = self.nts_score_sd
                      self.RNA["reactivityVector"] = self.reactivityVector
                      self.RNA["isENANtsV"] = ncm_by_nt_d_Dict["isENANtsV"]
                      
                      sc_so = [0]*len(ncm_by_nt_d_Dict["so"][0])
                      for ntTabBySS in ncm_by_nt_d_Dict["so"]:
                        for i in range(0,len(ntTabBySS)):
                          sc_so[i] += ntTabBySS[i]["scoreConfiance"]
                          
                      sc_mcff = [0]*len(ncm_by_nt_d_Dict["mcff"][0])
                      for ntTabBySS in ncm_by_nt_d_Dict["mcff"]:
                        for i in range(0,len(ntTabBySS)):
                          sc_mcff[i] += ntTabBySS[i]["scoreConfiance"]  
                      
                      self.RNA["sc_so"] = sc_so
                      self.RNA["sc_mcff"] = sc_mcff
                      
                      if(options.verbose):
                        print("############################ Apres le traitement des ARN et le calcul de leur score") 
                      stop7 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop7 - start))
                      #print("stn : "+str(e["stn"]))
                      if not os.path.exists(self.publicFolder):
                        os.makedirs(self.publicFolder)
                      
                      if(options.verbose):
                        print("file: "+self.publicFolder+str(self.id_ARN)+".json")
                      if(not self.error):
                        f = open(self.publicFolder+str(self.root)+"_"+str(self.id_ARN)+".json",'w')
                        f.write(json.dumps(self.RNA))
                        f.close()
                      else:
                        print("file error : "+self.filePath+str(self.id_ARN))
                      if(options.verbose):
                        print("############################ FIN $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$44") 
                      stop8 = timeit.default_timer()
                      if(options.verbose):
                        print (str(stop8 - start))
                  else:
                      print("no SS ... So: "+str(len(self.sESTabsubOpt))+ " MCFF :" +str(len(self.sESTabMcff)))

          self.fileError.close()
         
      def defVar(self,d,linesParsed):
        
          self.error = False
          self.d = d
          self.testNCM = 4
          
          self.filePath  = d["filePath"]
          self.expName = d["currentExp"]
          
          self.name = ""+str(self.expName)
          self.seqProfilTab = linesParsed
          self.RNA = {}
          self.RNA.clear()
          self.RNA_Str = ''
          
          
          self.indexTabAdded = []
          self.mcffTreshold = d["mcffTreshold"]
          self.soTreshold = d["soTreshold"]
        
          #pour donner une id unique a chaque sous-optimaux (structure secondaire)
          self.indexTabSO = {}
          self.indexTabMcff = {}
          
          #p = pca.PCA( array, fraction=fraction )
      
          self.id_ARN = str(d["uniqueId"])
          self.publicFolder = d["publicFolder"]
          self.prediction = d["prediction"]
          
          self.adenineCutOff = d["adenineCutOff"]
          self.scoreCutOff = d["scoreCutOff"]
          
          self.collectionName = d["collectionName"]
          self.dbName = d["dbName"]
          self.root = d["root"]
          self.MFP = d["MFP"]
      

        
      def toToWL(self,d):
        return " ".join(d.keys())

      def nT_Treatment(self,seq,ncm_by_nt_d_Dict,e,frP_mcff,frP_so,nt_id_num,mfp_mcff,mfp_so):
        nTs = []
        
        if("so" in ncm_by_nt_d_Dict and len(ncm_by_nt_d_Dict["so"][0]) > 0 and "mcff" in ncm_by_nt_d_Dict  and len(ncm_by_nt_d_Dict["mcff"][0]) > 0):
          for i in range(0,len(seq)):
            #print("ncm_by_nt_d_Dict" + str(ncm_by_nt_d_Dict["so"][i]))
            ncm_by_nt_d_Dict["so_u"] = [[]]*len(ncm_by_nt_d_Dict["so"])
            ncm_by_nt_d_Dict["mcff_u"] = [[]]*len(ncm_by_nt_d_Dict["mcff"])
            if(i >= len(ncm_by_nt_d_Dict["mcff_u"])):
              print("i :"+str(i)+ "len(ncm_by_nt_d_Dict['mcff_u'] : "+str(len(ncm_by_nt_d_Dict["mcff_u"])))
              self.error = True
              self.fileError.write("i >= len(ncm_by_nt_d_Dict_mcff,"+str(self.root)+";"+str(self.id_ARN)+"\n")
              return nTs
            if(i >= len(ncm_by_nt_d_Dict["so_u"])):
              print("i :"+str(i)+ "len(ncm_by_nt_d_Dict['so_u'] : "+str(len(ncm_by_nt_d_Dict["so_u"])))
              self.error = True
              self.fileError.write("i >= len(ncm_by_nt_d_Dict_so,"+str(self.root)+";"+str(self.id_ARN)+"\n")
              return nTs
            ncm_by_nt_d_Dict["so_u"][i] = [ncm_by_nt_d_Dict["so"][i][0]]
            ncm_by_nt_d_Dict["mcff_u"][i] = [ncm_by_nt_d_Dict["mcff"][i][0]]
            voisinPairAllsub_so = [x["voisins"] for x in ncm_by_nt_d_Dict["so_u"][i]]
            voisinPairAllsub_mcff = [x["voisins"] for x in ncm_by_nt_d_Dict["mcff_u"][i]]
            #print("voisinPairAllsub_so[0]  : "+str(voisinPairAllsub_so[0]))
            #print("voisinPairAllsub_mcff[0] : "+str(voisinPairAllsub_mcff[0]))
            nt_id_num = nt_id_num + 1
            score = e["scoreTab"][i]
            erreur = e["errorTab"][i]
            sorte = e["seq"][i]
            
            if(not(i< 1 or i > len(seq)-2)):
              seqMotif3 = e["seq"][i-1:i+2]
            else:
              seqMotif3 = "-"
            if(not(i< 2 or i > len(seq)-3)):
              seqMotif5 = e["seq"][i-2:i+3]
            else:
              seqMotif5 =  "-"
            
            freqPairee_mcff = frP_mcff[i]
            freqPairee_so = frP_so[i]
            position = i
            id_str = str(nt_id_num)
            reactivity_pred = self.reactivityVector[i]
            #
            ncm_Pred_bySS_so = ncm_by_nt_d_Dict["so_u"][i]
            ncm_Pred_bySS_mcff = ncm_by_nt_d_Dict["mcff_u"][i]
            inHeliceFinal_mcff_temp = 0
            inHeliceFinal_so_temp = 0
            
            if(i>1):
                #savoir si le nucleotide est au milieu d'une suite de 3 nucleotide pairee
                if(nTs[i-2]["freqPairee_so"] == 1 and freqPairee_so == 1 ):
                  #print("Hier than :"+ str(InH_sens))
                  inHeliceFinal_so = 1
                else:
                  inHeliceFinal_so = 0
                  
                if(nTs[i-2]["freqPairee_mcff"] == 1 and freqPairee_mcff == 1 ):
                  #print("Hier than :"+ str(InH_sens))
                  inHeliceFinal_mcff = 1
                else:
                  inHeliceFinal_mcff = 0
              
              
                voisinG = self.createVoisin(nTs[i-1])
                nt_i = i
                l = len(seq)
                region = self.findRegion(nt_i,l)
                nt = self.createNt(id_str,
                                    voisinG,
                                    erreur,
                                    score,
                                    sorte,
                                    position,
                                    self.id_ARN,
                                    inHeliceFinal_so_temp,
                                    inHeliceFinal_mcff_temp,
                                    freqPairee_mcff,
                                    freqPairee_so,
                                    seqMotif3,
                                    seqMotif5,
                                    mfp_mcff[i],
                                    mfp_so[i],
                                    region,
                                    voisinPairAllsub_so,
                                    voisinPairAllsub_mcff,
                                    [x["predColor"] for x in ncm_by_nt_d_Dict["so"][i]],
                                    [x["predColor"] for x in ncm_by_nt_d_Dict["mcff"][i]],
                                    reactivity_pred,
                                    ncm_Pred_bySS_so,
                                    ncm_Pred_bySS_mcff,
                                    [x["scoreConfiance"] for x in ncm_by_nt_d_Dict["so"][i]],
                                    [x["scoreConfiance"] for x in ncm_by_nt_d_Dict["mcff"][i]],
                                    ncm_by_nt_d_Dict["isENANtsV"])
                
                voisinD = self.createVoisin(nt)
                nTs[i-1]["voisinDroit"] = voisinD
                
                
                  
                nTs[i-1]["inHelice_so"] = inHeliceFinal_so
                nTs[i-1]["inHelice_mcff"] = inHeliceFinal_mcff
            else:
                voisinG = None
                inHeliceTemp = 0
            
            l = len(seq)
            region = self.findRegion(i,l)   
            #print("ncm_by_nt_d_Dict : "+str([x["predColor"] for x in ncm_by_nt_d_Dict["so"][i]]))
            nt = self.createNt(id_str,
                                voisinG,
                                erreur,
                                score,
                                sorte,
                                position,
                                self.id_ARN,
                                inHeliceFinal_so_temp,
                                inHeliceFinal_mcff_temp,
                                freqPairee_mcff,
                                freqPairee_so,
                                seqMotif3,
                                seqMotif5,
                                mfp_mcff[i],
                                mfp_so[i],
                                region,
                                voisinPairAllsub_so,
                                voisinPairAllsub_mcff,
                                [x["predColor"] for x in ncm_by_nt_d_Dict["so"][i]],
                                [x["predColor"] for x in ncm_by_nt_d_Dict["mcff"][i]],
                                reactivity_pred,
                                ncm_Pred_bySS_so,
                                ncm_Pred_bySS_mcff,
                                [x["scoreConfiance"] for x in ncm_by_nt_d_Dict["so"][i]],
                                [x["scoreConfiance"] for x in ncm_by_nt_d_Dict["mcff"][i]],
                                ncm_by_nt_d_Dict["isENANtsV"])
            nTs.append(nt)
            #si le nucleotide est le dernier, il n y a pas de voisin droit
            if (i == len(seq)-1):
                nTs[i]["voisinDroit"] = None
          
          #nTs_Hi = []
          #nTs_Low = []
          #nTs_Bg = []
          #for nt in nTs:
              ##print("nt : "+str(nt))
              #if(nt["score"] > self.d['hiScore']):
                  #nTs_Hi.append(nt)
              #else: 
                  #if(nt["score"] < self.d["lowScore"] 
                                          #and nt["score"] != -999):
                      #nTs_Low.append(nt)
                  #else:
                      #if(nt["score"] != -999):
                          #nTs_Bg.append(nt)
          
        return nTs
            
            
      #fin du traitement des nucleotides individuel -----------

      def createNt(self,
                   id_str,
                   voisinG,
                   erreur,
                   score,
                   sorte,
                   position,
                   rna_id,
                   inHelice_so,
                   inHelice_mcff,
                   freqPairee_mcff,
                   freqPairee_so,
                   seqMotif3,
                   seqMotif5,
                   mfp_mcff,
                   mfp_so,
                   region,
                   voisinPairAllsub_so,
                   voisinPairAllsub_mcff,
                   prediction_so,
                   prediction_mcff,
                   reactivity_pred,
                   ncm_Pred_bySS_so,
                   ncm_Pred_bySS_mcff,
                   sc_mcff,
                   sc_so,
                   ss_v):
          globalRNAinfo = {}
          d_temp = {}
          d_temp["u_id"] = id_str,
          d_temp["voisinGauche"] = voisinG
          d_temp["erreur"] = erreur
          d_temp["score"] = score
          d_temp["sorte"] = sorte
          d_temp["position"] = position
          d_temp["prediction_color_so"] = prediction_so
          d_temp["prediction_color_mcff"] = prediction_mcff
          d_temp["rna_id"] = rna_id
          d_temp["inHelice_so"] = inHelice_so
          d_temp["inHelice_mcff"] = inHelice_mcff
          d_temp["freqPairee_mcff"] = freqPairee_mcff
          d_temp["freqPairee_so"] = freqPairee_so
          d_temp["seqMotif3"] = seqMotif3
          d_temp["seqMotif5"] = seqMotif5
          d_temp["mfp_mcff"] = mfp_mcff
          d_temp["mfp_so"] = mfp_so
          d_temp["region"] = region
          d_temp["voisinPairAllsub_so"] = voisinPairAllsub_so
          d_temp["voisinPairAllsub_mcff"] = voisinPairAllsub_mcff
          d_temp["reactivity_pred"] = reactivity_pred
          d_temp["ncmTabDG_so"] = ncm_Pred_bySS_so
          d_temp["ncmTabDG_mcff"] = ncm_Pred_bySS_mcff
          d_temp["sc_mcff"] = sc_mcff
          d_temp["sc_so"] = sc_so
          d_temp["isENANts"] = ss_v[position]
          
          if(len(sc_mcff) > 0):
            d_temp["sc_mean_mcff"] = sum(sc_mcff)/len(sc_mcff)
          else:
            d_temp["sc_mean_mcff"] = 0 
            
          if(len(sc_so) > 0):
            d_temp["sc_mean_so"] = sum(sc_so)/len(sc_so)
          else:
            d_temp["sc_mean_so"] = 0  
          
          return d_temp
          
     
      
      
      
      def findMostfreqPartner(self,tabPairTab,nbSS):
          t1 = []
          if(len(self.seq) > 0):
            for i in range(0,len(self.seq)):
              t1.append([0]*len(self.seq))
            for pT in tabPairTab:
              for p in pT:
                #print("p : "+str(p))
                #print("len : "+str(len(tabPairTab)))
                if(len(p)>1):
                  if(p[1] > len(self.seq)):
                    if(self.options.verbose):
                      print("p[0] : "+str(p[0])+" p[1] : "+str(p[1]))
                      print("self.seq : "+str(len(self.seq)))
                      print("tabPairTab : "+str(tabPairTab))
                  t1[p[0]][p[1]] += 1
                  t1[p[1]][p[0]] += 1
                else:
                  break
            
            for i in range(0,len(self.seq)):
              for j in range(0,len(self.seq)):
                t1[i][j] = t1[i][j]/nbSS
              
          t2 = []
          for nt in t1:
            pairP = self.compileFreq(nt)
            #print("pairP : "+str(pairP))
            t2.append(pairP)
          
          
          return (t1,t2) 
              
      def compileFreq(self,tab):
        maxi = 0
        newK = -1
        for i in range(0,len(tab)):
          #print("tab[i] : "+str(tab[i]))
          if(tab[i] > maxi):
            newK = i
            maxi = tab[i]
        if(maxi == 0 ):
          newK = -1
          
        return newK
        
      def convertDBTo01(self,c):
        if(c == '.'):
          return 0
        else:
          return 1
        
      def sumNtWise(self,tabOfStruct):
        t = []
        if(len(tabOfStruct) > 0):
          t = [0]*len(tabOfStruct[0])
          for s in tabOfStruct:
            for i in range(0,len(s)):
              t[i] += self.convertDBTo01(s[i])
          for i in range(0,len(t)):
            t[i] = t[i]/len(tabOfStruct)
        return t
      
      
      def filterMinus999(self,tab):
        newTab = []
        for x in tab:
          if(x != -999):
            newTab.append(x)
        return newTab
      

      def findRegion(self,nt_i,l):
          if(nt_i/l < 0.05):
            return "_5p"
          if(nt_i/l > 0.95):
            return "_3p"
          return "millieu"
      
          
      def createVoisin(self,voisin):
          d_temp = {}
          d_temp["u_id"] = voisin["u_id"]
          d_temp["erreur"] = voisin["erreur"]
          d_temp["score"] = voisin["score"]
          d_temp["sorte"] = voisin["sorte"]
          return d_temp
      
      
      def calculFreqR(self,nt,seq):
          seq = seq.upper()
          freqR = seq.count(nt) / float(len(seq))
          return freqR
      
      def extractMotif(self,n,seq):
          alphabets = ['A', 'C', 'U', 'G']
          keywords = itertools.product(alphabets, repeat = n)
          tab = ["".join(x) for x in keywords ]
          d_temp = {}
          for e in tab:
              d_temp[e] = seq.count(e)
          return d_temp
      
      
      def jsNodeAll(self,tab,n):
          return [self.jsNode(pair[0],pair[1],1) for pair in tab]+self.phosphoLink(n)
      
      def jsNode(self,source,target,value):
          d = {}
          d["source"] = source
          d["target"] = target
          d["value"] = value
          return d

      def phosphoLink(self,n):
          tab = [{}]*(n-1)
          for i in range(0,n-1):
                  tab[i] = {}
                  tab[i]["source"] = i
                  tab[i]["target"] = i+1
                  tab[i]["value"] = 2
          return tab

      #Parser // Mcff (MCFlashfold) RNAso / RNAfold (Vienna package) need to be instaled on your system

      def mcff(self,seq,adjMat):
          
          if (1):#(len(seq) < self.mcffTreshold):
                  lineSplited = []
                  self.isThereAnAux = 1
                  stri = 'mcff'+str(" -s " + "'"+seq.upper()+ "' -ft "+str(self.mcffTreshold)+ " -tables /u/malricp/mcff/MC-Flashfold-v34/src/")
                  if(self.options.verbose):
                    print ("mcffresult : "+stri)
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  lineSplited = output.strip().decode('ascii').split()
                  #print("len lineSplited mcff:"+str(len(lineSplited)))
                  self.sESTabMcff = self.parseMcff(lineSplited)
                  #structureEnergieTab = structureEnergieTab[:80]
                  if(1):#(len(seq) < self.mcffTreshold):
                      self.adjMatMcff = self.adMat2(self.sESTabMcff)
                  self.adjListMcff = self.createAdList(self.indexTabMcff,self.adjMatMcff)

          return self.sESTabMcff


      def parseMcff(self,lineSplited):
          structureEnergyTab = []
          struct = ""
          rank = 0
          #self.sESTabMcff = []
          
          for index in range(0,len(lineSplited)):
                  if(lineSplited[index] == "bad"):
                      structureEnergyTab = []
                      break
                  #construction du tableau de sous-optimaaux
                  i=index%3
                  if i == 0:
                      #
                      struct = lineSplited[index]
                      if(self.options.verbose):
                        print ("lineSplited[index] : "+str(lineSplited[index]))
                        print("struct : "+struct)
                        print("structLen : "+str(len(struct)))
                      #if ('(' not in struct):
                              #break
                  elif i== 1 :
                      energy = lineSplited[index]
                  elif i == 2:
                      #print(energy)
                      summary = lineSplited[index]
                      sE = StructEnergy(struct,self.seq,float(energy),summary,rank)
                      #print "strucEnergy McFF : " + sE.displayStructEnergy()
                      structureEnergyTab.append(sE)
                      self.indexTabMcff[struct] = math.floor(index/3)
                      rank = rank+1
          return structureEnergyTab[:self.mcffTreshold]

      def rnaSO(self,seq,adjMat):
          if (1):#len(seq)<self.soTreshold:
                  
                  lineSplited = []
                  structureEnergieTab = []
                  #r=10
                  #self.currentrange = (r*self.step,r+self.numStep*self.step)
                  #print "range: ["+str(r*self.step)+','+str(self.smallestTab[r])+']'
                  self.stableRTab = []
                  self.stableCutedTab = []
                  if not os.path.exists(self.filePath+"fa/"):
                          os.makedirs(self.filePath+"fa/")
                  if not os.path.exists(self.filePath+"fa/"):
                          os.makedirs(self.filePath+"fa/")
                  fastaP = self.filePath +"fa/"+ self.name+".fa"
                  #print ("fastaP = "+fastaP)
                  self.fastaWrite(fastaP,seq)
                  #self.namew = self.name+str(r)
                  stri = "RNAsubopt -e "+str(self.d["so_e_value"])+" -s < \""+fastaP+"\""
                  #print (stri)
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  #print ("out : "+output.strip().decode('ascii'))
                  lineSplited = output.strip().decode('ascii').split()
                  #la sequence
                  #print lineSplited[3]
                  #print ("line de SO : "+str(len(lineSplited)))
                  if(len(lineSplited)>3):
                          self.seq = lineSplited[3]
                          structureEnergieTab = self.parseSO(lineSplited[6:])
                          #print "ParseSo"
                          #structureEnergieTab = structureEnergieTab[:40]
                          #print("self.sESTabsubOpt: "+ "|".join([str(x) for x in self.sESTabsubOpt]))
                          #print("structureEnergieTab len : "+str(len(structureEnergieTab)))
                          if(1):#(len(seq) < self.mcffTreshold):
                              self.adjMatSO = self.adMat2(structureEnergieTab)
                          self.createAdjListSO(self.filePath+"AdjList")
                          #printPair(self.pairTabSo)

                          #----------------------
                          #GraphMaker_compare(self.seqTab,self.pairTabTab,self.name)#-----------------------

                          #print ("stable region : "+str(self.stableR))
                  return structureEnergieTab
        

      def parseSO(self,lineSplited):
          finalTab = []
          rnaSuboptTab = []
          for i in range(0,len(lineSplited)):
                      #print("i : "+str(i)+" line : "+str(lineSplited[i]))
                      if ( i%2 == 1):
                          #print("1 : " + lineSplited[i])
                          e = lineSplited[i]
                          rnaSuboptTab.append(st+" "+str(e))
                      else:
                          #print("0 : " + lineSplited[i])
                          st = lineSplited[i]
          #print ("nombre de SO : "+str(len(rnaSuboptTab)))
          if len(rnaSuboptTab)>0:
                      #print "mfe_SO : "+rnaSuboptTab[0]
                      #print "#SO : "+str(len(rnaSuboptTab))
                      
                      for i in range(0,len(rnaSuboptTab)):
                          #print ("rnaSuboptTab["+str(i)+"] : "+rnaSuboptTab[i])
                          rsoSplited = rnaSuboptTab[i].split()
                          if (not is_number(rsoSplited[1]) and float(rsoSplited[1]) == 0):
                                  #print("breakOn : "+rsoSplited[0])
                                  self.error = True
                                  self.fileError.write(str(self.root)+";"+str(self.id_ARN)+"\n")
                                  break
                          sE = StructEnergy(rsoSplited[0],self.seq,float(rsoSplited[1]),"-",i)
                          #print("####rsoSplited[0] : "+rsoSplited[0])
                          self.indexTabSO[rsoSplited[0]] = i
                          #print ("strucEnergy SO : " + sE.displayStructEnergy())
                          finalTab.append(sE)
          return finalTab[:self.soTreshold]
        

       # I use a file (fasta) to input RNAso and RNAfold


      def fastaWrite(self, outFile,seq):
        if not os.path.exists(outFile+"fa/"):
            try:  
                    with open(outFile,"w") as newFasta:
                            header = ">" + self.name + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                            newFasta.write(header + seq)
            except IOError:
                    print ("Failed to open " + outFile)
                    exit(1)

        #print ">> Done!"



      def findPairGen(self,secStruct):
          openDict = {}
          openDict['('] = 1 
          openDict['['] = 1 
          openDict['{'] = 1 
          openDict['<'] = 1 
          openDict['A'] = 1 
          openDict['B'] = 1 
          tempTab=[]
          #print self.SEnergyTabMcff[0].struct
          stList = list(secStruct)
          #Pour tous les nucleotides (, trouve sa paire
          for n in range(0,len(stList)):
                  #print("recherche de la pair du nucleotide:"+str(n))
                  c = stList[n]
                  if c in openDict :
                          i = self.findPartnerGen(n,secStruct)
                          if(i == -1):
                                  print ("findPAirNOMA : " + str(n))
                          tempTab.append((n,i))
          return tempTab
          
          
          
          

      def pbToNcmAll(self,detailed,mfe):
        # trouve le NCM de toutes les SS et predit le score de reactivite. 
        # detailled est une variable booleen qui est vrai si on veut les MCN detaillezz
        # detaillez = avec la sequence et la position
        l = len(self.seq)
        ntTab_so = [[] for x in range(0,l)]
        if(mfe):
          l_so = 1
        else:
          l_so = len(self.pairTabAllSubOpt_so)
        ncmSStab = []
        for j in range(0,l_so): #pour tous les sous-optimaux (so)
          oneSS_so = self.pbToNcm(self.pairTabAllSubOpt_so[j],l,"so",detailed)
          ncmSStab.append(oneSS_so["ncmTabD"])
          for i in range(0,l):
            #print("i : "+str(i)+" len tab :"+str(len(oneSS_so["ncmPredictionNtTab"])))
            nt = {}
            nt.clear()
            if(len(oneSS_so["ncmPredictionNtTab"]) > 0):
              nt["predColor"] = oneSS_so["ncmPredictionNtTab"][i]
              nt["scoreConfiance"] = oneSS_so["scoreTab"][i]
            nt["ncm"] = oneSS_so["ncmTabD"][i]
            nt["voisins"] = oneSS_so["voisinTab"][i]
            ntTab_so[i].append(nt)
            
        isEachNt_AllNcmthesameVector = self.eachNt_AllNcmthesameVector(ncmSStab)
            
        
        ntTab_mcff = [[] for x in range(0,l)]
        if(mfe):
          l_mcff = 1
        else:
          l_mcff = len(self.pairTabAllSubOpt_mcff)
        ntTab_mcff = [[] for x in range(0,l)]
        for j in range(0,l_mcff):#pour tous les sous-optimaux (mcff)    
          oneSS_mcff = self.pbToNcm(self.pairTabAllSubOpt_mcff[j],l,"mcff",detailed)
          for i in range(0,l):
            nt = {}
            nt.clear()
            if(len(oneSS_mcff["ncmPredictionNtTab"]) > 0):
              nt["predColor"] = oneSS_mcff["ncmPredictionNtTab"][i]
              nt["scoreConfiance"] = oneSS_mcff["scoreTab"][i]
            nt["ncm"] = oneSS_mcff["ncmTabD"][i]
            nt["voisins"] = oneSS_mcff["voisinTab"][i]
            ntTab_mcff[i].append(nt)
        
        
        return_d = {}
        return_d["isENANtsV"] = isEachNt_AllNcmthesameVector
        return_d["so"] = ntTab_so #[x[0] for x in ntTab_so]
        return_d["mcff"] = ntTab_mcff #[x[0] for x in ntTab_mcff]
        
        return return_d
      
      def eachNt_AllNcmthesameVector(self,ncmSStab):
        firstTab = [x["merge"] for x in ncmSStab[0]]
        returnTab = [1]*len(ncmSStab[0])
        missmatchTab = [0]*len(ncmSStab[0])
        for ss_i in range(0,len(ncmSStab)):
          first = ncmSStab[ss_i][0]["merge"]
          if(self.options.verbose):
            print("first : "+first)
          missmatch = 0
          #pour chaque nucleotide de chaque Structure secondaire
          for nt_i in range(1,len(ncmSStab[ss_i])):
            if(firstTab[nt_i] != ncmSStab[ss_i][nt_i]["merge"]):
              missmatchTab[nt_i] +=1
              if(missmatchTab[nt_i] > 0):
                #print("ncmSStab[ss_i][nt_i]['merge'] : "+ ncmSStab[ss_i][nt_i]["merge"] +" != "+firstTab[nt_i])
                returnTab[nt_i] = 0
        #print("returnTab : "+",".join(str(x) for x in returnTab))
        return returnTab
      
      
      # fonction qui transforme un tableau de paire de bases en MCN
      def pbToNcm(self,paireBaseTab,l_seq,soft,detailed):
              
              #print(",".join("("+str(x[0])+","+str(x[1])+")" for x in paireBaseTab))
              if(self.prediction):
                self.db = self.openMongoClient()
              #Le dictionnaire des pair de base (d[origine] = cible) se tableau fonctionne des 2 bords
              d_pb = self.getPB(paireBaseTab)
              ncmTabD = []
              
              #print("paireBaseTab: "+";".join("("+str(x[0])+"|"+str(x[1])+")" for x in paireBaseTab))

              voisinTab = l_seq*[None]
              for i in range(0,l_seq):
                vd = self.voisinDroit(i,l_seq,d_pb)
                vg = self.voisinGauche(i,d_pb)
                voisinTab[i] = self.getVoisin(d_pb,i,l_seq,vd,vg)
              
              
              #pour tous les nucleotides
              for i in range(0,l_seq):
                  vd = self.voisinDroit(i,l_seq,d_pb)
                  vg = self.voisinGauche(i,d_pb)
                  #print("s : "+s)
                  ncmTabD.append({})
                  ncmTabD[i].clear()
                   
                  if(vd >= 0 and vg >= 0 and vd == int(d_pb[vg])):
                      #literallement : (le premier voisin de droite pairee) egale (le nt pairee avec le premier voisin de gauche pairee)
                      #print("--res :"+str(len(range(vd,int(d_pb[vd])+1))))
                      s = str(str(vd-d_pb[vd]+1)+self.seqNcm("-",self.seq[vg:vd+1],detailed)+self.posNcm(str(i-vg),True))
                      ncmTabD[i]['merge'] = s
                      #loop
                  else:
                      if(i in d_pb):
                          #nucleotide Paire
                          ncmTabD[i] = self.treatNtPaired(i,d_pb,vd,vg,detailed)
                                  
                      else:
                          #le ncm n'est pas pairee
                          self.treatNtNonPaired(ncmTabD,d_pb,vd,vg,i,detailed)    
                                
                  if('droit' not in ncmTabD[i]):
                      ncmTabD[i]['droit'] = "-"
                      
                  if('gauche' not in ncmTabD[i]):
                      ncmTabD[i]['gauche'] = "-"
                      
                  if('merge' not in ncmTabD[i]):
                      ncmTabD[i]['merge'] = ncmTabD[i]['gauche']+"&"+ncmTabD[i]['droit']
                  
              
              ncmPredictionNtTab = []
              scoreTab = []
              #print("ncmTabD : "+str(len(ncmTabD)))
              for j in range(0,len(ncmTabD)):
                if(self.prediction):
                  (color,scoh) = self.getPredColorweigth(ncmTabD[j]["merge"],j,soft)
                  #if(win == 'red' or win == 'blue'):
                    #print("color : "+win)
                  ncmPredictionNtTab.append(color)
                  scoreTab.append(scoh)
                else:
                  ncmPredictionNtTab.append("black")
                  scoreTab.append(1)
              if(self.prediction):
                self.db.close
              ntDict = {}
              ntDict["voisinTab"] = voisinTab
              #For all nt we have a dict with the key : droit, gauche et merged
              ntDict["ncmTabD"] = ncmTabD
              # the color of the reactivity prediction (Blue = hi, red = Low, grey = no prediction)
              ntDict["ncmPredictionNtTab"] = ncmPredictionNtTab
              # the sum of the scoreTab (divide by the number of prediction) is a measure of how well the structure is acurate
              ntDict["scoreTab"] = scoreTab
              return ntDict


   

      def seqNcm(self,seq1,seq2,detailed):
        if(detailed):
          if(len(seq1) > 6):
            seq1 = "L"
          if(len(seq2) > 6):
            seq2 = "L"
          return "-"+seq1+"-"+seq2
        else:
          return ""
      
      def posNcm(self,pos,detailed):
        if(detailed):
          return "_pos_"+str(pos)
        else:
          return ""
      
      
      def openMongoClient(self):
        d = {}
        client = MongoClient()
        client = MongoClient('10.0.0.161', 27027)
        db = client[self.dbName]
        return db
      
      def getPrediction2(self,low,bg,hi):
        l_low = low
        l_bg = bg
        l_hi = hi
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_bg / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Bg"
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_hi / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Hi"
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_low / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Low"
        
        return "NED"
      
      def getPrediction(self,low,bg,hi):
        l_low = low
        l_bg = bg
        l_hi = hi
        if((l_hi + l_low)>0):
          poids = abs(l_hi - l_low)/(l_hi + l_low)
          if((l_low + l_hi ) > self.d['nbNt_Ncm'] and l_hi / float(l_low + l_hi ) > self.d["predCutOff"] ):
            return ("Hi",poids)
          if((l_low + l_hi ) > self.d['nbNt_Ncm'] and l_low / float(l_low + l_hi ) > self.d["predCutOff"] ):
            return ("Low",poids)
        return ("NED",0)
      
      def colorTranlation(self,pred):
        d = {"Hi":"blue","Low":"red","NED":"white","Bg":"purple"}
        return d[pred]
      
      def getColor2Ncm(self,p1,p2):
        if(p1 ==  p2):
          return self.colorTranlation(p1)
        if(p1 == "Bg" or p1 == "NED"):
          return self.colorTranlation(p2)
        if(p2 == "Bg" or p2 == "NED"):
          return self.colorTranlation(p1)
        # p1 and p2 are Hi and Low or Low and Hi ...
        return "indigo"
      
      # Donne la couleur du contour des MCN et cacul le score de prediction
      def getPredColorweigth(self,ncm,i,soft):
        if(len(self.reactivityVector) > i):
            actualLevel = self.reactivityVector[i]
            scoh = 0
            color = "black"
            if(actualLevel == "Bg"):
              color = "grey"
            else:
              #print("ncm : "+ ncm)
              cursor = self.db[self.collectionName].find({"ncm":ncm,"soft":soft})
              isNcmpresent = False
              for d in cursor:
                isNcmpresent = True
                low = d["low"]
                #print("low : "+ str(low))
                hi = d["hi"]
                #print("hi : "+ str(hi))
                bg = d["bg"]
                #print("bg : "+ str(bg)+"\n")
              if(isNcmpresent):
                (pred,poids) = self.getPrediction(low,bg,hi)
                color = self.colorTranlation(pred)
              else :
                color = "green"

            if(color == "red" and actualLevel == "Low" or color == "blue" and actualLevel == "Hi"):
              scoh = 1 * poids
            if(color == "blue" and actualLevel == "Low" or color == "red" and actualLevel == "Hi"):
              scoh = -1 * poids

            return (color,scoh)
        else:
            self.error = True
            self.fileError.write("reactivityVector too small,"+str(self.root)+";"+str(self.id_ARN)+"\n")
            return ("yellow", 0)
      
      
      def getPB(self,paireBaseTab):
        d = {}
        for i in range(0,len(paireBaseTab)):
            #print('value :'+str(paireBaseTab[i]))
            t = paireBaseTab[i][0]
            s = paireBaseTab[i][1]
            d[s] = t
            d[t] = s
            #print('s :' + str(s) + ' t = '+str(t))
        return d
      
      
      def getVoisin(self,d_pb,i,l_seq,vd,vg):
        
        d = {}
        d['i'] = i 
        
        d['vd'] = vd
        if(vd in d_pb):
          d['p_vd'] = d_pb[vd]
        else:
          d['p_vd'] = "-"
        
        
        d['vg'] = vg
        if(vg in d_pb):
          d['p_vg'] = d_pb[vg]
        else:
          d['p_vg'] = "-"
          
         
        if(vd >= 0 and vg >= 0 and vd == int(d_pb[vg])):
            if(i not in d_pb):
              d['p_i'] = -1
            else:
              d['p_i'] = d_pb[i]
        
        if(i not in d_pb):
            d['p_i'] = -1
        else:
            d['p_i'] = d_pb[i]
            
        return d
      
      
      
      
      def voisinDroit(self,i,l,d_pb):
          if(i < l):
              #print("vg i : "+str(i))
              for j in range(i+1,l):
                  #print("j : "+str(j))
                  if(j in d_pb):
                      return j
              return -1
          else:
              return -2

      def voisinGauche(self,i,d_pb):
            if(i > 0):
                #print("vd i : "+str(i))
                for j in range(i-1,-1,-1):
                    #print("j : "+str(j))
                    if(j in d_pb):
                        return j
                return -1
            else:
                return -2
  
      def treatNtPaired_droit(self,i,d_pb,vd,detailed):
          if(i == d_pb[vd]):
              loop_v = vd-i+1
              if(loop_v>10):
                  s = "L"
              else:
                  s = str(str(loop_v)+self.seqNcm("-",self.seq[i:vd+1],detailed)+self.posNcm(str(0),True))
                  
              #loop
          else:
              ntBetweVd1 = vd-i+1
              ntBetweVd2 = d_pb[i]-d_pb[vd]+1
              #print("vg1 : "+ str(vd)+" vg2 : "+ str(d_pb[vd])+" i1 : "+ str(i)+" i2 : "+ str(d_pb[i]))
              #print("--res :"+str(ntBetweVd1)+','+str(ntBetweVd2))
              if(abs(ntBetweVd1) > 6 ):
                  ntBetweVd1Str = "L"
                  seqBt1a = "L"
                  seqBt1b = "L"
              else:
                  ntBetweVd1Str = str(ntBetweVd1)
                  seqBt1a = ntBetweVd1Str
                  seqBt1b = self.seq[i:vd+1]
                  
              if(abs(ntBetweVd2) > 6 ):
                  ntBetweVd2Str = "L"
                  seqBt2a = "L"
                  seqBt2b = "L"
              else:
                  ntBetweVd2Str = str(ntBetweVd2) 
                  seqBt2a = ntBetweVd2Str
                  seqBt2b = self.seq[d_pb[vd]:d_pb[i]+1]
                  
              s = seqBt1a+'_'+seqBt2a+self.seqNcm(seqBt1b,seqBt2b,detailed)+self.posNcm(str(0),True)
              
              
          return s
      
      
      
      def treatNtPaired_gauche(self,i,d_pb,vg,detailed):
          
          if(i == d_pb[vg]):
              loop_v = (i-vg)+1
              if(loop_v>10):
                s = "L"
              else:
                s = str(str(loop_v)+self.seqNcm("-",self.seq[vg:i+1],detailed)+self.posNcm(str(i-vg),True))
              
              #loop
          else:
              ntBetweVg1 = i-vg+1
              ntBetweVg2 = d_pb[vg]-d_pb[i]+1
              #print("vg1 : "+ str(vg)+" vg2 : "+ str(d_pb[vg])+" i1 : "+ str(i)+" i2 : "+ str(d_pb[i]))
              #print("--res :"+str(ntBetweVg1)+','+str(ntBetweVg2))
              if(ntBetweVg1 > 6 or -1*ntBetweVg1 > 6):
                  ntBetweVg1Str = "L"
              else:
                  ntBetweVg1Str = str(ntBetweVg1)
                  
              if(ntBetweVg2 > 6 or -1*ntBetweVg2 > 6):
                  ntBetweVg2Str = "L"
              else:
                  ntBetweVg2Str = str(ntBetweVg2) 
                  
              s = str(ntBetweVg1Str+'_'+ntBetweVg2Str + self.seqNcm(self.seq[vg:i+1],self.seq[d_pb[i]:d_pb[vg]+1],detailed)+self.posNcm(str(i-vg),True))
              
              
          return s
      
      
      
      def treatNtPaired(self,i,d_pb,vd,vg,detailed):
          d = {}
          if(vd >= 0 ):
              #VoisinDroitPaire
              d['droit'] = self.treatNtPaired_droit(i,d_pb,vd,detailed)
              
          else:
              if(vd == -1):
                  d['droit'] = "noVD_Paired"
              if(vd == -2):
                  d['droit'] = "Dernier_Paired"
                  
          if(vg >= 0 ):
              d['gauche'] = self.treatNtPaired_gauche(i,d_pb,vg,detailed)
          else:
              if(vg == -1):
                  d['gauche'] = "noVG_Paired"
              if(vg == -2):
                  d['gauche'] = "Premier_Paired"
                  
          return d
      
      
      def treatNtNonPaired(self,ncmD,d_pb,vd,vg,i,detailed) :
        if(vd < 0 or vg < 0):
            if(vd == -1):
                ncmD[i]['droit'] = "noVD_notPaired"
            if(vg == -1):
                ncmD[i]['gauche'] = "noVG_notPaired"
            if(vd == -2):
                ncmD[i]['droit'] = "Dernier_notPaired"
            if(vg == -2):
                ncmD[i]['gauche'] = "Premier_notPaired"
        else:
            ntBetweVdEtVg1 = vd-vg+1
            ntBetweVdEtVg2 = d_pb[vg]-d_pb[vd]+1
            #print("vg1 : "+ str(vg)+" vg2 : "+ str(d_pb[vg])+" vd1 : "+ str(vd)+" vd2 : "+ str(d_pb[vd]))
            #print("--res :"+str(ntBetweVdEtVg1)+','+str(ntBetweVdEtVg2))
            if(ntBetweVdEtVg1 > 6 or -1*ntBetweVdEtVg1 > 6):
                ntBetweVdEtVg1Str = "L"
            else:
                ntBetweVdEtVg1Str = str(ntBetweVdEtVg1)
                
            if(ntBetweVdEtVg2 > 6 or -1*ntBetweVdEtVg2 > 6):
                ntBetweVdEtVg2Str = "L"
            else:
                ntBetweVdEtVg2Str = str(ntBetweVdEtVg2)    
            
            s = str(ntBetweVdEtVg1Str+'_'+ntBetweVdEtVg2Str + self.seqNcm(self.seq[vg:vd+1],self.seq[d_pb[vd]:d_pb[vg]+1],detailed)+self.posNcm(str(i-vg),True))
            ncmD[i]['nonPaired'] = s
            ncmD[i]['merge'] = s
      
      



      def addD(self,d1,d2):
          dFin = {}
          for (key,value) in d1.items():
              if key in d2 :
                  dFin[key] = value + d2[key]
              else:
                  dFin[key] = value
          for (key,value) in d2.items():
              if not (key in d1) :
                  dFin[key] = value
          return dFin
          

      def findPartnerGen(self,n,secStruct):
          stList = list(secStruct)
          oppose2 = createDictOpenToClose()
          oppose = createDictCloseToOpen()
          acc = {}
          for key in oppose2:
                  acc[key] = 0
          #En realit/ on a besoin que d'un compteur
          
          stList = list(secStruct)
          acc[stList[n]] = 1
          for i in range(n+1,len(stList)):
                  if stList[i] != '.':
                          if stList[i] == stList[n]:
                                  acc[stList[i]] += 1
                          if stList[i] in oppose :
                                  if stList[n] == oppose[stList[i]]:
                                              acc[oppose[stList[i]]] -= 1
                                              #print ("index : "+str(i)+" -> "+str(acc))
                          if acc[stList[n]] == 0 :
                                  return i
          print ("no match2 : " +stList[n])
          print (secStruct + " : "+str(n))
          return -1
        

        
      def output_files (self,t):
          self.namew = self.name.split('/')[0]
          print ("NAME *****************: " +self.jsFileout+t+".js")
          print ("stat/"+self.specificPath)
          with open(self.jsFileout+t+".js", 'w') as f :
                  f.write(self.jsOut)
                  f.close
          with open(self.htmlFileout+t+".html", 'w') as f :
                  f.write(self.htmlOut)
                  f.close

      def findMotif(self,motif):    
          mo = [m.start() for m in re.finditer(motif,self.seq)]
          iMotif = []
          for m in mo:
                  for k in range(0,len(motif)):
                          iMotif.append(m+k)
          return iMotif


      def jsDictGeneratorByPredictor(self,mcffStr,soStr):
          d_temp = {}
          d_temp["mcff"] = mcffStr
          d_temp["rnaSubOpt"] = soStr
          return d_temp
          


      def nodeT(self,nom,energy,predictor,p_int):
          d = {}
          d["name"] = nom
          d["predictor"] = predictor
          d["energy"] = energy
          d["p_int"] = p_int
          return d


      def linkT(self,i1,i2,value,pred):
          d = {}
          d["source"] = i1
          d["target"] = i2
          d["value"] = value
          d["pred"] = pred
          return d

      def charTogroup(self,c):
          c = c.upper()
          #print str(c)
          group = ''
          if c == 'A':
                  group = 0
          if c == 'T':
                  group = 1
          if c == 'U':
                  group = 1
          if c == 'C':
                  group = 2
          if c == 'G':
                  group = 3
          return group




  #-------------------------------------------fastaMulti


      def fastaMultiTos(self,filePath):
          inFile  = filePath
          outFile = inFile +".out" 
          if(options.verbose):
            print (">> Opening FASTA file...")
          # Reads sequence file list and stores it as a string object. Safely closes file:
          try:  
                  with open(inFile,"r") as newFile:
                          sequences = newFile.read()
                          sequences = re.compile("^>", flags=re.MULTILINE).split(sequences) # Only splits string at the start of a line.
                          del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
                                          # Del removes this empty element.
                          newFile.close()
          except IOError:
                  print ("Failed to open " + inFile)
                  exit(1)

          if(options.verbose):
            print (">> Converting FASTA file from multiline to single line and writing to file.")
          # Conversts multiline fasta to single line. Writes new fasta to file.
          for fasta in sequences:
                  try:
                          header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
                  except ValueError:
                          print ("fasta : "+fasta)
                  self.numberOf_windows = int(math.floor(len(sequence)/self.step))
                  for r in range (0,self.numberOf_windows-self.numStep):
                          with open(outFile+str(r)+".fa","w") as newFasta:
                                  if(options.verbose):
                                    print ("lenSeq : "+str(len(sequence))+" r : "+str(r))
                                  header_new = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                                  s = sequence.replace("\n","")
                                  sequence = s.replace('\r','') + "\n" # Replace newlines in sequence, remember to add one to the end.
                                  smallest = len(sequence)
                                  if ((r+self.numStep)*self.step < len(sequence)):
                                          smallest = ((r+self.numStep)*self.step)
                                  self.smallestTab.append(smallest)
                                  if(options.verbose):
                                    print ("r to smallest : ("+str(smallest)+","+str((r+self.numStep)*self.step)+")")
                                  newFasta.write(header_new + sequence[(r*self.step):smallest])
                                  newFasta.close()
                                  if(options.verbose):
                                    print (">> Done! "+outFile+str(r)+".fa")
        
      #processing



      def adMat2(self,sESTab):
          # 'd' contiendra un dictionnaire de dictionnaire 
          d = {}
          #pour toute les sous-optimaux, on fait la liste des structure voisine (exponentiel ... sur la longueur de la sequence)
          #print("adMat skipped!")
          #sESTab = []
          for so in sESTab:
                      soPairTab1 = self.findPairGen( so.struct )
                      lPairTab = len(soPairTab1)
                      if(lPairTab >0):
                              d[so.struct] = {}
                              for so2 in sESTab:
                                  soPairTab2 = self.findPairGen( so2.struct )
                                  lPairTab = max(len(soPairTab2),lPairTab)
                                  ressemblance = self.commonPair(soPairTab1,soPairTab2)/float(lPairTab)
                                  #print ("resemblance = "+str(ressemblance))
                                  if(ressemblance > 0.95):
                                          #if(so.struct in self.indexTabMcff and so2.struct in self.indexTabMcff):
                                              #print (str(self.indexTabMcff[so.struct])+","+str(so2.energy - so.energy)+","+str(self.indexTabMcff[so2.struct]))
                                          d[so.struct][so2.struct]     =    ressemblance
                        #print "Total : "+so.struct +" : "+ str(len(dtemp)) + " sur : "+ str(len(sESTab))
          return d

      def commonPair(self,so1Tab,so2Tab):
          count = 0
          for p1 in so1Tab:
                      for p2 in so2Tab:
                              if(p1[0] == p2[0] and p1[1] == p2[1]):
                                  #print ("common pair :"+str(p1[0])+","+str(p2[0])+","+str(p1[1])+","+str(p2[1]))
                                  count += 1
          return count

      def createAdList(self,indexTab,adjMat):
          tab = []
          for s in adjMat:
                      #print "-----"+name
                      for s2 in adjMat[s]:
                              if(s != s2):
                                #print (str(str(indexTab[s])+";"+str(adjMat[s][s2])+";"+str(indexTab[s2])))
                                if (-1 *adjMat[s][s2]< 0):
                                    #f.write(s+","+str(adjMat[s][s2])+","+s2+"\n")
                                    tab.append((indexTab[s],adjMat[s][s2],indexTab[s2]))
          #f.close()
          return tab

        
      def createAdjListSO(self,path):
          self.adjListSO = []
          if not os.path.exists(path+"/adList/"):
                      os.makedirs(path+"/adList/")
         
          for s in self.adjMatSO:
                      #print "-----"
                      for s2 in self.adjMatSO[s]:
                              if(s != s2):
                                first = str(self.indexTabSO[s])
                                second = str(self.adjMatSO[s][s2])
                                third = str(self.indexTabSO[s2])
                                #print ("createAdjListSO : "+first+";"+second+";"+third)
                                if (-1 * self.adjMatSO[s][s2]<0):
                                    if(s in self.indexTabSO and s2 in self.indexTabSO):
                                            #f.write(s+","+str(self.adjMatSO[s][s2])+","+s2+"\n")
                                            self.adjListSO.append((self.indexTabSO[s],self.adjMatSO[s][s2],self.indexTabSO[s2]))
       
    
      #mAdToD3Graph prend en entree deux listes d'adjencences et deux listes de conformations et donne un graph de transition en format noeud + link
      def mAdToD3Graph(self,sEtabMCff,sEtabSO,adjListMcff,adjListSO):
          maxMcff = max([x.energy for x in sEtabMCff])
          minMcff = min([x.energy for x in sEtabMCff])
          maxSO = max([x.energy for x in sEtabSO])
          minSO = min([x.energy for x in sEtabSO])
          nodeTempsT = ""
          tabNode = []
          for i in range(0,len(sEtabMCff)):
                      tabNode.append(self.nodeT(i,sEtabMCff[i].energy,"mcff",1))

          for i in range(0,len(sEtabSO)):
                      tabNode.append(self.nodeT(i+len(sEtabMCff),sEtabSO[i].energy,"so",2))

          linkTempsT = ""
          tabLink = []
          if(self.options.verbose):
            print("adjList len : "+str(len(sEtabMCff)))
          for triplet in adjListMcff:
                      tabLink.append( self.linkT(triplet[0],triplet[2],triplet[1],"mcff"))#(source,target,value)
                              
          for triplet in adjListSO:
                      tabLink.append( self.linkT(triplet[0]+len(sEtabMCff),triplet[2]+len(sEtabMCff),triplet[1],"so"))#(source,target,value)

          for k in self.sESTabMcff:
                      s1 = k.struct
                      for k2 in  self.sESTabsubOpt:
                                s2 = k2.struct
                                #print("s1:"+s1+"\ns2:"+s2)
                                soPairTab2 = self.findPairGen( s2 )
                                soPairTab1 = self.findPairGen( s1 )
                                lPairTab = max(len(soPairTab2),len(soPairTab1))
                                ressemblance = self.commonPair(soPairTab1,soPairTab2)/float(lPairTab)
                                if(ressemblance > 0.8):
                                          print("2 structures semblable!")
                                          #print("so : "+str(self.indexTabSO[s2]))
                                          #print("mcff : "+str(self.indexTabMcff[s1]))
                                          tabLink.append( self.linkT(self.indexTabSO[s2]+len(sEtabMCff),self.indexTabMcff[s1],1,"so_mcff"))#(source,target,value)

          
          #print ("tabLink : "+str(tabLink[:100]))
          return (tabNode,tabLink,maxMcff,minMcff,maxSO,minSO)


      def displayCount(self):
          print ("Nombre de structure sous-optimales de MCff: %d et de RNAsubOpt : %d"  % (len(self.SEnergyTabMcff),len(self.sESTabsubOpt)))

      def restartCount():
          RNA2Daio.count=0

      def displayRNA2D(self):
          print ("Sequence : ", self.seq)
          if len(self.SEnergyTab) != 0  :
                  for s in self.SEnergyTab:
                          #print "type =",type(s)
                          s.displayStructEnergy()
          else :
                  print ("longueur Energy tab nul")
      

    
def normalise(d,kind):
    if(kind == "sum"):
      s = sum(d.values())
      if(s == 0):
          #print("sum = 0")
          #for i in d.itervalues():
            #print(str(i))
          factor = 0
      else:
          factor=1.0/s
      for k in d:
          d[k] = d[k]*factor
    return d

def createDictOpenToClose():
    oppose2 = {}
    oppose2['(']=')'
    oppose2['{']='}'
    oppose2['[']=']'
    oppose2['<']='>'
    oppose2['A']='a'
    oppose2['B']='b'
    return oppose2 

def createDictCloseToOpen():
    oppose = {}
    oppose[')']='('
    oppose['}']='{'
    oppose[']']='['
    oppose['>']='<'
    oppose['a']='A'
    oppose['b']='B'
    return oppose



def printPair(ptab):
  for p in ptab:
          print ("("+str(p[0])+","+str(p[1])+")")





#-----------------------------------------------------------------------------------------------------------------------------------------------StructureEnergyTuple
class StructEnergy:
    'structure/Energy'
    count = 0

    def __init__(self, struct,seq, energy, summary,rank):
          self.seq = seq
          self.struct = struct
 
          self.energy = energy
          self.summary =  summary
          self.rank =  rank
          StructEnergy.count += 1
          
          
    def __str__(self):
          return "\nSeq  : "+self.seq+"\nSSec : "+self.struct + "\nEnergy : "+ self.energy
          
    def toString(self):
          return self.struct + " "+ str(self.energy)
    

    def displayCount(self):
          print ("Nombre de structure 2D %d" % StructEnergy.count)
  
    def restartCount(self):
          StructEnergy.count=0

    def displayStructEnergy(self):
          s = "".join(["Structure: ", self.struct,  ", Energy: ", str(self.energy), " Summary: ", self.summary, " Rank: ", str(self.rank)])
          return s
  
  #processing


    def findPartnerFromStruct(self,n):
          acc = 1
          stList = list(self.struct)
          for i in range(n+1,len(stList)):
                  if stList[i] == '(':
                          acc += 1
                  if stList[i] == ')':
                          acc -= 1
                          #print ("index : "+str(i)+" -> "+str(acc))
                  if acc == 0 :
                          return i
          print ("no match")
          print (self.struct + " : "+str(n))
          return -1







iparens = iter('()')
parens = dict(zip(iparens, iparens))
closing = parens.values()

