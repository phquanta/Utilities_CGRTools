from typing import List
from selenium import webdriver
import time

from numpy.core.fromnumeric import prod, sort
import rdkit
from rdkit.Chem import AllChem
from rdkit import Chem, RDConfig
from rdkit.Chem import MolFromSmiles, AllChem
from rdkit.Chem.rdChemReactions import PreprocessReaction
import CGRtools as cgr
from pathlib import Path
import numpy as np
from io import StringIO
import os
import argparse
from icecream import ic
from CGRtools.reactor import *

from SmilesEnumerator import SmilesEnumerator
import itertools
import collections
from random import shuffle

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
import subprocess
import xml.etree.ElementTree as ET

#def getSMILESCGRS(*args):


from checkEnumerator import CGREnumerator

class formSMIRKSromLst(object):
        def __init__(self,inRx,outRx):
            self.inRx= inRx
            self.outRx=outRx
        def __str__(self):
            return ".".join(self.inRx)+">>"+".".join(self.outRx)
        def __repr__(self):
                return ".".join(self.inRx)+">>"+".".join(self.outRx) 

    
class SeleniumMarvinJs(object):
    def __init__(self,timing1=3,timing2=4):
        self.timing1=timing1
        self.timing2=timing2
        self.browser = webdriver.Firefox()
        self.browser.get('http://mapper.grzybowskigroup.pl/marvinjs/')
        self.textarea = self.browser.find_element_by_id("wejscieUzytkownika")
        self.btnMapReaction=self.browser.find_element_by_xpath("//input[@type='submit' and @value='Map reaction!']")
        self.btnMapDrawing=self.browser.find_element_by_id("zmapujRysunek")
        self.iframe=self.browser.find_element_by_xpath("//iframe")

    def getProperReactionSMILES(self,smiles)->str:
        self.textarea.clear()
        self.textarea.send_keys(smiles)
        self.btnMapReaction.click()
        time.sleep(self.timing1)
        self.textarea.clear()
        self.browser.switch_to.frame(self.iframe)
        self.browser.find_element_by_css_selector("[title='Add/Remove explicit H']").click()
        self.browser.switch_to.default_content()
        time.sleep(self.timing1)
        self.btnMapDrawing.click()
        time.sleep(self.timing2)
        properSMILES=self.textarea.get_property('value')
        return properSMILES
        
        
        
    

class CGRUitls(object):
    
    
    @staticmethod
    def dumpListToFile(lst,file)->None: 
         with open(file, 'w') as f:
            for i,rx in enumerate(lst):
                f.write(f"{rx}\n")
                
    @staticmethod
    def readFileToList(file)->List: 
        lst=[]
        with open(file, 'r') as f:
            for line in f: 
                 line = line.strip() 
                 lst.append(line) 
        return list(set(lst))
                
        
        
    def __init__(self,useSelenium=True):
        self.sme = SmilesEnumerator()
        if useSelenium:
            self.selenium=SeleniumMarvinJs()
        else:
            self.selenium=None
        self.enumerator=CGREnumerator()   
        #enumerator.enhanceCGR(rx)
        return
    
    
    def getSMIRKS(self,file)->None: 
        self.lstSMIRKS=self.readFileToList(file)
        print(self.lstSMIRKS)
    def canonicalize(self,lst)->list:
        rx=[]
        for i,r in  enumerate(lst):
           container=next(cgr.files.SMILESRead(StringIO(r)))
           #ic(type(container))
           container.canonicalize()
           rx.append(str(container))
        return rx
    
    def createVariation(self,lst)->list:
        shuffle(lst)
        return [self.sme.randomize_smiles(rx) for rx in lst]
        #return [rx for rx in lst]

    def getCGRsGivenSMIRKS(self,rxLst)->list:
        cgrs=[]
        if not os.path.exists("./results/worker.log"):
            try:
                os.mkdir("./results")
                with open("./results/worker.log", 'a'):
                    os.utime("./results/worker.log", None)
            except OSError:
                print("Error in Creating worker.log")
            
            
            
        
            
        for i,rx in enumerate(rxLst):
            ic("###############################################")
            ic(rx)
            if True:
               process=subprocess.Popen(["map_reaction",f"{rx}"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
               stdout,stderr = process.communicate()
            
               root = ET.fromstring(stdout)
               rxAAM=""
               for el in root.findall("./solutions/solution/reactants"):
                    rxAAM=el.text
                    
               for el in root.findall("./solutions/solution/products"):
                        rxAAM=rxAAM+">>"+el.text

#            ic(type(rx))
#            if len(rx)==0: 
 #               continue
            #rxAAM=str(rx)
                        
#                     ic(child.tag, child.attrib)
            #ic(str(rx))
            #rxAAM=rxAAM.replace(")-[",")[")
            #rxAAM=rxAAM.replace(")-(",")(")
            #rxAAM=rxAAM.replace("]-[","][")
            
            #rxAAM=rxAAM.replace("(-[","([")
            
            #s1 = AllChem.ReactionFromSmarts(rxAAM)
            #rxAAM=AllChem.ReactionToSmarts(s1)
            
            
            #rxAAMmarvin=rxAAM
            rxAAMarvin=rxAAM
            for _ in range(1):
                rxAAMmarvin=self.selenium.getProperReactionSMILES(rxAAMarvin)
                
                
                ic(rxAAM)
                ic(rxAAMmarvin)
                if len(rxAAMmarvin)==0:
                        continue
                try:
                
                    reactionObj=next(cgr.files.SMILESRead(StringIO(rxAAMmarvin)))
#                    print("Length REACTIONOBJ:"))
                    #reactionObj.clean2d(randomize=True)
                    ic(type(reactionObj))
                    reactionObj.standardize()
                    reactionObj.implicify_hydrogens()
                    reactionObj.clean2d()
                    ic("SEEEEEEEEEEEee",str(reactionObj))
                    
                except Exception as e:
                    print("Error, yet continuing")    
                    print(e)
                
                
                
                cgrObj=~reactionObj    
                print(str(cgrObj))
                #if True:
                if len(cgrObj.centers_list)>0:
                    ic(str(reactionObj),type(reactionObj))

                

                    #cgrObj.clean2d()
                    ic(str(cgrObj))
                    ic(len(cgrObj.centers_list))
                    ic(cgrObj.centers_list)
                #    if len(cgrObj.centers_list)>1:
                #        continue
                    if(str(cgrObj) not in cgrs):
                        cgrs.append(str(cgrObj))
                        
                    ic("Length-CGR,Rx Number:",len(cgrs),i)    
            
            ic(len(rxLst))
            ic(len(cgrs))

        return cgrs

            
            
            
        
    
    def createRxVariatons(self,*args)->list:
        inRx,outRx,numVariations=args
        rxLst=[]
        for i in range(numVariations):
            rxLst.append(formSMIRKSromLst(self.createVariation(inRx),self.createVariation(outRx)))
        rxLst=list(set(rxLst))
        return rxLst    
               
            
                
        
        
        
    
    def enhanceCGRs(self,numVariations,rx)->list:
        rxAll=[]
        self.enumerator.nTrials=numVariations
        self.enumerator.enhanceCGR(rx)
        self.enumerator.unwind=False
        rxAll=self.enumerator.cgrS[:]
        
#          for i in enumerator.cgrS:
 #               ic(i)
        
        return rxAll
            
         
    
    def isValidReaction(self,reactionObj)->None:
        valid=False
        reactionMols=list(reactionObj.reactants)
        reactionMols.extend(reactionObj.products)
        reactionMols.extend(reactionObj.reagents)
        
        for i in reactionMols:
            try:
                i.kekule()
                
                
                lst=i.check_valence()
                if len(lst)==0 and i.check_thiele():
                    valid=True
                
            except Exception as e:
                ic("ERROR IN isValidReactrion")
                return False
        return valid   


    def standartize(self,rx)->None:
            rx.standardize()
            rx.kekule()
            rx.implicify_hydrogens()
            rx.clean2d()



    def cgr2SMIRKS(self,cgrFile)->[list,int,int]:

        generator=cgr.files.SMILESRead(cgrFile)
        reactionsSMIRKS=set()
        cntValid=0  
        cntAll=0 
        while True:
            try:
                    nextReactionObj=next(generator)
                    if type(nextReactionObj)==cgr.containers.cgr.CGRContainer: 
                        cntAll+=1
                        cgrContainer=nextReactionObj
                        rx = cgr.ReactionContainer.from_cgr(cgrContainer)
                        self.standartize(rx)
                        
                    elif type(nextReactionObj)==cgr.containers.ReactionContainer:    
                        rx=nextReactionObj
                        self.standartize(rx)

                    if self.isValidReaction(rx):
                                cntValid+=1
                                reactionsSMIRKS.add (str(rx))
                            
            except StopIteration:
                ic("########### STOP ITERATION ##################")
                break
            
            except Exception as e:
                ic("########### SOME ERROR ##################")
                
                print(e)

        return [list(reactionsSMIRKS),cntValid,cntAll]
