

from subprocess import Popen,PIPE



def parseMcff(lineSplited):
    structureEnergyTab = []
    rank = 0
    #self.sESTabMcff = []
    
    for index in range(0,len(lineSplited)):
            if(lineSplited[index] == "bad"):
                structureEnergyTab = []
                break
            #construction du tableau de sous-optimaaux
            i=index%3
            if i == 0:
                #print ("lineSplited[index] : "+str(lineSplited[index]))
                struct = lineSplited[index]
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
                rank = rank+1
    return structureEnergyTab[:10]
  
  
  
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
  



stri = "mcff -s 'GGGCGGGAAACGGAUUGGAUUAAGAAGGAGGCACUUAAUCCGAUCCGCUAAGCGGCGCAAGCCUCCAGCGCCGUAUGUCGCCCAAGCGGUGCUUCGGCACCGCAAAAGAAACAACAACAACAAC' -ft 1"
print (stri)
Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)

output,err =  Process.communicate()

print("output : "+output.strip().decode('ascii'))

print("err : "+err.strip().decode('ascii'))


lineSplited = output.strip().decode('ascii').split()
print("len lineSplited mcff:"+str(len(lineSplited)))
sESTabMcff = parseMcff(lineSplited)


