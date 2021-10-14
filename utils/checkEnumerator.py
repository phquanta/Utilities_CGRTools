from logging import exception
from builtins import Exception, dict, object
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
import re
from CGRtools.containers import *
import CGRtools as cgr
from io import StringIO
import numpy as np
from numpy import asarray
import random    





class CGREnumerator(object):
    _instance = None
    _MAX_FAILURES=500
    _MAX_SUCCESSIVE=100
    _OFFSET_TOCOUNT_MEAN=10
    _DEFAULT_NUM_TRIALS=1
    _ATOMS = sorted([
                'Al', 'As', 'B', 'Br', 'C', 'Cl',  'K', 'Li', 'N',
                'Na', 'O', 'P', 'S', 'Se', 'Si', 'Te','Mg','Cu','Zn','Pd','Pt','Fe',
                'Pb','Hg','Zr','Mn','Ti','Ca','Cr','Sb','Ce','Rh','Ru',
                'Ag','Tl','Re','Ir','Ni','Ba','Au','Ge','W','Mo','Xe','Ga','Cd','Ta','Bi','He',
                'Sm','V','Nd','Be','Sr', 'Dy','Yb','La','Y'
            ],key=len)



    _CGR_CORE_GRAMMAR=['[.>-]', '[.>=]', '[.>#]', '[.>:]', '[.>~]', '[->.]', '[->=]', '[->#]', '[->:]',\
                '[->~]', '[=>.]', '[=>-]',  '[=>#]', '[=>:]', '[=>~]', '[#>.]', '[#>-]', '[#>=]',\
                '[#>:]', '[#>~]', '[:>.]', '[:>-]', '[:>=]', '[:>#]',  '[:>~]', '[~>.]', '[~>-]',\
                '[~>=]', '[~>#]', '[~>:]']
    _CGR_PROBLEM_SPECIFIC=['[N+>0]','[O->0]','[O0>-]','[Fe+>+3]','[N0>+]','Cl->0])',\
                 '[I->0]','[Mn+>+3]']
    _NOT_IMPLEMENTED=[]
    
    #_NOT_IMPLEMENTED=[ '[N+>0]','[O->0]','[N0>+]','[O0>-]']
    
    
    _CGRGRAMMAR=_CGR_CORE_GRAMMAR+_CGR_PROBLEM_SPECIFIC

    _REGEXPRESSION='(?<!\[|-|=|\)|\>)\.'
    
    """
    def __new__(cls):
        if cls._instance is None:
            print('Creating the object')
            cls._instance = super(CGREnumerator, cls).__new__(cls)
            cls._instance.nTrials=cls._DEFAULT_NUM_TRIALS
            cls._instance.cgrS=[] 
            cls._instance.sme= SmilesEnumerator()  
            cls._instance._rxInverse=[]
            cls._instance.successRatio=0.
            cls._instance.eps=1.E-5
            cls._instance.unwind=False
        return cls._instance
    """    
    
    def __init__(self,nTrials=1)->None:
            self.nTrials=self._DEFAULT_NUM_TRIALS
            self.sme= SmilesEnumerator()  
            self.eps=1.E-5
            self.unwind=False
        
     

    


    def createVariation(self,lst)->list:
        
        shuffle(lst)
        try:
            return [self.sme.randomize_smiles(rx) if random.choice([True,False])==True else rx for rx in lst ]
        except Exception as e:
            return [rx  for rx in lst ]
    
    def enhanceCGR(self,rx,Flag=True)->None:
        self.cgrS=[]
        self._rxInverse=[]
        counter=0
        #for _ in range(self.nTrials):
        counterFailed=0
        lengths=[]
        
        while True:    
            counter=counter+1
            
            lstSMILES=[]
            compounds=re.split(self._REGEXPRESSION,rx)
            
            #ic(compounds)
            
            compounds_copy=compounds[:]
            for elem in compounds_copy:
                if not any([x in elem for x in self._CGRGRAMMAR]):
                    lstSMILES.append(elem)
                    compounds.remove(elem)
        
            compounds_copy=compounds[:]   
            compounds=[]
            for elem in compounds_copy:
                    cgrIdx=np.where(asarray([x in elem for x in self._CGRGRAMMAR])==True)[0]
                    notCGrIdx=np.where(asarray([x in elem for x in self._CGRGRAMMAR])==False)[0]
                    cgrBonds=list(asarray(self._CGRGRAMMAR)[cgrIdx])
                    cgrBond_original=cgrBonds[:]
                    
                    
                    
                    idx=[random.randrange(0,len(notCGrIdx)) for x in range(len(cgrIdx))]
                    atomsPresent=asarray(self._ATOMS)[np.where(asarray([x in elem for x in self._ATOMS])==True)[0]]
                    atomsNotPresent=asarray(self._ATOMS) \
                                         [np.where(asarray([x in elem for x in self._ATOMS])==False)[0]][0:len(cgrIdx)]
                    atomsNotPresent=["[*:"+str(i)+"]" for (i,x) in enumerate(atomsNotPresent)]
                    dictReplacement=dict(zip(cgrBonds,atomsNotPresent))
                    rx_new=elem
                    if not any([x in rx_new for x in self._NOT_IMPLEMENTED]) and Flag:
#                        if Flag:
                            for key, value in dictReplacement.items():
                                
                                rx_new=rx_new.replace(key,value)
                                
                            #ic(elem,rx_new)
                            
                            try:
                                rx_new=self.sme.randomize_smiles(rx_new)
                            except Exception as e:
                                continue    
                            for key, value in dictReplacement.items():
                                rx_new=rx_new.replace(value,key)

                    compounds.append(rx_new)
                    
            del compounds_copy    

            _order=random.choice([True,False])
            shuffle(compounds)
            rxSMILES=self.createVariation(lstSMILES)
        
                
            cgrSTR=".".join(1 and (rxSMILES.extend(compounds) if _order else compounds.extend(rxSMILES) )\
                                or (rxSMILES if _order else compounds))
            
            
            try:
                cgrObj=next(cgr.files.SMILESRead(StringIO(cgrSTR)))
                self.decomposed = ReactionContainer.from_cgr(cgrObj)
                #ic(str(self.decomposed))
                if cgrSTR not in self.cgrS:
                    self.cgrS.append(cgrSTR)
                    #ic(cgrSTR)
            except Exception as e:
                #ic("2",e)
                  #ic(counterFailed,Flag,counter)
                  counterFailed=counterFailed+1
                  if counterFailed>self._MAX_FAILURES and Flag==True:
                        self.unwind=True
                        self.enhanceCGR(rx,Flag=False) 
                  elif Flag==False and  counterFailed>=self._MAX_FAILURES:  
                        return
                  else:
                        continue  
            if self.unwind==True and Flag==True: return        
#            ic("HERE",counterFailed,Flag,counter,self.cgrS)
            self.decomposed.clean2d()
            
            if str(self.decomposed) not in self._rxInverse:
                self._rxInverse.append(str(self.decomposed))
            
            self.successRatio=len(self.cgrS)/counter
            if len(self.cgrS)==self.nTrials: return
            #ic(counterFailed)
            if self.successRatio<self.eps and Flag==False and  counterFailed>self._MAX_FAILURES:
                return
            lengths.append(len(self.cgrS))
            if len(lengths)>self._OFFSET_TOCOUNT_MEAN:
                #ic(lengths[-self._OFFSET_TOCOUNT_MEAN:])
                #ic(np.mean(lengths[-self._OFFSET_TOCOUNT_MEAN:]))
                if int(np.mean(lengths[-self._OFFSET_TOCOUNT_MEAN:]))==len(self.cgrS) and counter>=self._MAX_SUCCESSIVE:
                    return
            
if __name__ == '__main__':
        #rx="[N+]CCC(Al[N+]CCCC[N+])BrO"
        #rx="O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O"

        #rx="[P-](F)(F)(F)(F)(F)F.C([.>-]Nc1ccc(N2CCOCC2)nc1)(Cc3ccc(c4cc(ncc4)C)cc3)([->=]O)[=>.]O.c56c(cccn5)nnn6OC(N(C)C)=[N+](C)C.N(C=O)(C)C.N(CC)(C(C)C)C(C)C"
        #rx="[N+](Cc1ccccc1)(CCCC)(CCCC)CCCC.c2(C([.>-]C([->.][Na])#N)[->.]Cl)c(cc(cc2)CC(NC(c3c(N4CCCCC4)cccc3)CCC)=O)OCC.C(Cl)Cl.O.[K+].[Cl-].[I-]"
        #rx="S(c1ccc(C)cc1)(O[->.]C([.>-]Oc2c(c3ccccc3)cccc2)CCCCCCl)(=O)=O.N(C=O)(C)C.C(=O)([O-])[O-].O.[K+].[K+]"
        #rx="C(CC[N+][Ta+]C[Bk-]O)C([O-])=O"
        #rx="C(C(N[->.]C(CNC)=O)CCC(N)=O)(=O)O.CN"

        #rx="O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O"
        #rx="O[=>-]O.[N+>0][->.]C(CC[N+])[.>=]O"
        #rx="O[=>-]O.CC(C)CCCC(C)CCCC(CCCC(C(C(=O)[O-])[->=]O)C)C"
        #rx="O[=>-]O.[N+>0][->.]C(CC[N+])[.>=]O"
        #rx="O[=>-]O.O=C(C(C(O)C([->=]O)CO)O)CO"
        #rx="O[=>-]O.O.N"
        #rx="O[=>-]O.NCCCC([->.]N)[.>=]O"
        
        #rx="O[->.]O.O[->.]O"
        
        rx="C[.>-]1([->.]C([.>=]O[->.]O)(=O)[->.]C[.>-]1(=O)[O-])CCC[N+]"
        enumerator=CGREnumerator()   
        enumerator.nTrials=10
        enumerator.enhanceCGR(rx)
        
        for i in enumerator.cgrS:
            ic(i)
        
        ic(enumerator.successRatio)        
        ic(len(enumerator._rxInverse))
        ic(len(enumerator.cgrS))









