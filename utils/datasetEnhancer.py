# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 10:22:29 2020

@author: Andrei
"""


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
import selenium
from utils import CGRUitls


def main(baseFile,cgrFile,enhanceFile,fromCGRFile=True):
    ll=[]
    if not fromCGRFile:
        cgr=CGRUitls()
        cgr.getSMIRKS(baseFile)
        lst=cgr.enhanceSMIRKS(1)
        cgr.dumpListToFile(lst,enhanceFile)
        #lst=["[O:1]=[O:2].[H][O:3][C:4]([H])([O:5][H])[c:6]1([H])[o:7][c:8]([H])([C:9]([H])=[O:10])[c:11]([H])[c:12]1[H]>>[H][O:1][O:2][H].[H][C:9](=[O:10])[c:8]1([H])[o:7][c:6]([H])([c:12]([H])[c:11]1[H])[C:4]([O-:3])=[O:5]"]
        cgr.dumpListToFile(cgr.getCGRsGivenSMIRKS(lst),cgrFile)
    else:
        cgr=CGRUitls(useSelenium=False)
        lst=cgr.readFileToList(cgrFile)
        lstAll=[]
        for i,elem in enumerate(lst):
            ic(elem)
            
            
            ic("######################################")
            lCGR=cgr.enhanceCGRs(40,elem) 
            ll.append(len(lCGR))
            if len(lCGR)==0:
                lCGR=[elem]
            for j in lCGR:
                ic(j)
            
            print(f"\n {i} \n")    
        
        
        
        ic(ll)
        ic(np.sum(ll))
        
            
        
            
            
    


if __name__ == "__main__":
    

#
#    parser = argparse.ArgumentParser()
#    parser.add_argument("infile", type=str,
#                    help="Input Initial Filename with reduced data")
#    parser.add_argument("outfile",
#                    help="Output Reactions  SMILES with ")
#    args = parser.parse_args()
    
#    main(args.infile, args.outfile)
    #main("../FineTuneRx_OO.smi", "../out.smi","../outCGR_10x.smi")
    main("../FineTuneRx_OO.smi", "../outCGR_1x.smi","../outCGR_10x.smi")
    
