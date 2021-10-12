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
    _ATOMS = sorted([
                'Al', 'As', 'B', 'Br', 'C', 'Cl',  'K', 'Li', 'N',
                'Na', 'O', 'P', 'S', 'Se', 'Si', 'Te','Mg','Cu','Zn','Pd','Pt','Fe',
                'Pb','Hg','Zr','Mn','Ti','Ca','Cr','Sb','Ce','Rh','Ru',
                'Ag','Tl','Re','Ir','Ni','Ba','Au','Ge','W','Mo','Xe','Ga','Cd','Ta','Bi','He',
                'Sm','V','Nd','Be','Sr', 'Dy','Yb','La','Y'
            ],key=len)



    _CGRGRAMMAR=['[.>-]', '[.>=]', '[.>#]', '[.>:]', '[.>~]', '[->.]', '[->=]', '[->#]', '[->:]',\
                '[->~]', '[=>.]', '[=>-]',  '[=>#]', '[=>:]', '[=>~]', '[#>.]', '[#>-]', '[#>=]',\
                '[#>:]', '[#>~]', '[:>.]', '[:>-]', '[:>=]', '[:>#]',  '[:>~]', '[~>.]', '[~>-]',\
                '[~>=]', '[~>#]', '[~>:]']

    _REGEXPRESSION='(?<!\[|-|=|\)|\>)\.'
    
    def __init__(self,nTrials)->None:
        self.nTrials=nTrials
        self.cgrS=[] 
        self.sme= SmilesEnumerator()  
        self._rxInverse=[]
        self.succesRation=0.

    


    def createVariation(self,lst)->list:
        
        shuffle(lst)
        return [self.sme.randomize_smiles(rx) for rx in lst]
    
    def enhanceCGR(self,rx)->None:
        for _ in range(self.nTrials):
            lstSMILES=[]
            compounds=re.split(self._REGEXPRESSION,rx)
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
                    for key, value in dictReplacement.items():
                        rx_new=rx_new.replace(key,value)
                    try:
                        rx_new=sme.randomize_smiles(rx_new)
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
                decomposed = ReactionContainer.from_cgr(cgrObj)
                if cgrSTR not in self.cgrS:
                    self.cgrS.append(cgrSTR)
            except Exception as e:
                continue
            
            decomposed.clean2d()
            
            if str(decomposed) not in self._rxInverse:
                self._rxInverse.append(str(decomposed))
            
            self.successRatio=len(self.cgrS)/self.nTrials
            #ic(len(cgrS),len(_rxInverse))
            #ic(success_ratio)
            #ic()
            
if __name__ == '__main__':
        #rx="[N+]CCC(Al[N+]CCCC[N+])BrO"
        #O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O

        #rx="[P-](F)(F)(F)(F)(F)F.C([.>-]Nc1ccc(N2CCOCC2)nc1)(Cc3ccc(c4cc(ncc4)C)cc3)([->=]O)[=>.]O.c56c(cccn5)nnn6OC(N(C)C)=[N+](C)C.N(C=O)(C)C.N(CC)(C(C)C)C(C)C"
        #rx="[N+](Cc1ccccc1)(CCCC)(CCCC)CCCC.c2(C([.>-]C([->.][Na])#N)[->.]Cl)c(cc(cc2)CC(NC(c3c(N4CCCCC4)cccc3)CCC)=O)OCC.C(Cl)Cl.O.[K+].[Cl-].[I-]"
        rx="S(c1ccc(C)cc1)(O[->.]C([.>-]Oc2c(c3ccccc3)cccc2)CCCCCCl)(=O)=O.N(C=O)(C)C.C(=O)([O-])[O-].O.[K+].[K+]"
        #rx="C(CC[N+][Ta+]C[Bk-]O)C([O-])=O"
        #rx="C(C(N[->.]C(CNC)=O)CCC(N)=O)(=O)O.CN"

        #rx="O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O"
        enumerator=CGREnumerator(1222)   
        enumerator.enhanceCGR(rx)


        
        for i in enumerator.cgrS:
            ic(i)
        
        ic(enumerator.successRatio)        
        ic(len(enumerator._rxInverse))
        ic(len(enumerator.cgrS))









