from logging import exception
from builtins import Exception, dict
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

sme = SmilesEnumerator()

#rx="[N+]CCC(Al[N+]CCCC[N+])BrO"
#O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O

#rx="[P-](F)(F)(F)(F)(F)F.C([.>-]Nc1ccc(N2CCOCC2)nc1)(Cc3ccc(c4cc(ncc4)C)cc3)([->=]O)[=>.]O.c56c(cccn5)nnn6OC(N(C)C)=[N+](C)C.N(C=O)(C)C.N(CC)(C(C)C)C(C)C"
#rx="[N+](Cc1ccccc1)(CCCC)(CCCC)CCCC.c2(C([.>-]C([->.][Na])#N)[->.]Cl)c(cc(cc2)CC(NC(c3c(N4CCCCC4)cccc3)CCC)=O)OCC.C(Cl)Cl.O.[K+].[Cl-].[I-]"
rx="S(c1ccc(C)cc1)(O[->.]C([.>-]Oc2c(c3ccccc3)cccc2)CCCCCCl)(=O)=O.N(C=O)(C)C.C(=O)([O-])[O-].O.[K+].[K+]"
#rx="C(CC[N+][Ta+]C[Bk-]O)C([O-])=O"
#rx="C(C(N[->.]C(CNC)=O)CCC(N)=O)(=O)O.CN"

#rx="O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O"

atoms = sorted([
            'Al', 'As', 'B', 'Br', 'C', 'Cl',  'K', 'Li', 'N',
            'Na', 'O', 'P', 'S', 'Se', 'Si', 'Te','Mg','Cu','Zn','Pd','Pt','Fe',
            'Pb','Hg','Zr','Mn','Ti','Ca','Cr','Sb','Ce','Rh','Ru',
            'Ag','Tl','Re','Ir','Ni','Ba','Au','Ge','W','Mo','Xe','Ga','Cd','Ta','Bi','He',
            'Sm','V','Nd','Be','Sr', 'Dy','Yb','La','Y'
        ],key=len)

ic(atoms)

cgrGrammar=['[.>-]', '[.>=]', '[.>#]', '[.>:]', '[.>~]', '[->.]', '[->=]', '[->#]', '[->:]',\
             '[->~]', '[=>.]', '[=>-]',  '[=>#]', '[=>:]', '[=>~]', '[#>.]', '[#>-]', '[#>=]',\
              '[#>:]', '[#>~]', '[:>.]', '[:>-]', '[:>=]', '[:>#]',  '[:>~]', '[~>.]', '[~>-]',\
               '[~>=]', '[~>#]', '[~>:]']

regExpression='(?<!\[|-|=|\)|\>)\.'

stack=[]

mapping={"Ta+": "->.",
         "Rf": ".>-",
         "Bk-":".>=",
         
    
}

def createVariation(lst)->list:
        global sme
        shuffle(lst)
        return [sme.randomize_smiles(rx) for rx in lst]





m = re.split(regExpression, rx)
ic(rx)
print(m)




_rxInverse=[]


nn=22
for _ in range(nn):
    lstSMILES=[]
    compounds=re.split(regExpression,rx)
    compounds_copy=compounds[:]
    for elem in compounds_copy:
        #ic(len(compounds))
        #print([x in elem for x in cgrGrammar])
        
        #ic(cgrIdx)
        
        
        if not any([x in elem for x in cgrGrammar]):
        #if len(cgrIdx)==0:
            lstSMILES.append(elem)
            compounds.remove(elem)
        
    compounds_copy=compounds[:]   
    #ic(compounds_copy)
    compounds=[]
    for elem in compounds_copy:
            cgrIdx=np.where(asarray([x in elem for x in cgrGrammar])==True)[0]
            notCGrIdx=np.where(asarray([x in elem for x in cgrGrammar])==False)[0]
            #ic(cgrIdx)
            cgrBonds=list(asarray(cgrGrammar)[cgrIdx])
            cgrBond_original=cgrBonds[:]
            #cgrBonds=[]
            """
            for el in cgrBond_original:
                #a=[m.start() for m in re.finditer(el, elem)]
                a=elem.find(el)
                ic(a)
                ic(el)
                ic(len(elem))
                #ic(elem[10])
                adds=""
                if elem[a-1] not in ["]",")"] and elem[a+len(el)] not in ["[","("]:
                    adds=elem[a-1]+el+elem[a+len(el)]
                    ic("adds1",adds)
                elif elem[a-1] in ["]",")"] and elem[a+len(el)] in ["[","("]: 
                    adds=el    
                    ic(adds)
                else:    
                    adds=el+elem[a+len(el)]    
                    ic("adds2",adds)
                    
                #adds=elem[a-1]+el+elem[a+len(el)] if elem[a-1] not in ["]",")"] else 
                print(adds)
             #   cgrBonds.append(adds)
             """
            import random    
            idx=[random.randrange(0,len(notCGrIdx)) for x in range(len(cgrIdx))]
            atomsPresent=asarray(atoms)[np.where(asarray([x in elem for x in atoms])==True)[0]]
            #aa=np.random.randint(len)
            atomsNotPresent=asarray(atoms)[np.where(asarray([x in elem for x in atoms])==False)[0]][0:len(cgrIdx)]
            #atomsNotPresent=asarray(atoms)[np.where(asarray([x in elem for x in atoms])==False)[0]][idx]
            #atomsNotPresent=["["+x+"-]" for x in atomsNotPresent]
            #atomsNotPresent=["[*:"+str(i)+"]" for (i,x) in enumerate(atomsNotPresent)]
            atomsNotPresent=["(*:"+str(i)+")" for (i,x) in enumerate(atomsNotPresent)]
            #atomsNotPresent=["[*]" for (i,x) in enumerate(atomsNotPresent)]
            #atomsNotPresent=["" for (i,x) in enumerate(atomsNotPresent)]
            dictReplacement=dict(zip(cgrBonds,atomsNotPresent))
            ic(dictReplacement)
            ic(atomsNotPresent)
            rx_new=elem
            for key, value in dictReplacement.items():
                rx_new=rx_new.replace(key,value)
                pass
            ic("---------------------------------------")
            ic(elem)
            ic(rx_new)
            try:
                rx_new=sme.randomize_smiles(rx_new)
                pass
            except Exception as e:
                ic(e)
                continue    
            ic(rx_new)
            for key, value in dictReplacement.items():
                  rx_new=rx_new.replace(value,key)
                  pass
            ic(rx_new)
            compounds.append(rx_new)
            #compounds.append(elem)
               
                
            

            #ic(atomsPresent)
            #ic(atomsNotPresent)
            #ic(cgrBonds)
            #ic(dict(zip(cgrBonds,atomsNotPresent)))
            
            #ic(asarray(atoms)[atomsPresent])
            
    del compounds_copy    
    
    import random 
    _order=random.choice([True,False])
    _order=True
    shuffle(compounds)
    rxSMILES=createVariation(lstSMILES)
    
    ic(rxSMILES)
    
    
    print("###############################################")
    
    cgrSTR=".".join(1 and (rxSMILES.extend(compounds) if _order else compounds.extend(rxSMILES) )\
                        or (rxSMILES if _order else compounds))
    ic(cgrSTR)
    #r="O[=>-]O.C(CC[N+][->.]C[.>=]O)C([O-])=O"
    #container=next(cgr.files.SMILESRead(StringIO(r)))
    
    try:
        cgrObj=next(cgr.files.SMILESRead(StringIO(cgrSTR)))
        decomposed = ReactionContainer.from_cgr(cgrObj)
    except Exception as e:
        ic(e)     
        continue
    
    
    decomposed.clean2d()
    ic(str(decomposed))
    
    if str(decomposed) not in _rxInverse:
        _rxInverse.append(str(decomposed))
    ic(_rxInverse)
    ic(len(_rxInverse))
    
    ic(rx)
    ic(lstSMILES)        
    ic(compounds)        
    
    
            
        
    
"""
    rx_new=sme.randomize_smiles(rx)
    for key, value in mapping.items():
        rx_new=rx_new.replace(key,value)

    
    
    print(rx_new)
"""
