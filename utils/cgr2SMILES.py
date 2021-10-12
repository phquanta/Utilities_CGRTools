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
from utils import CGRUitls


def main(cgrFile,reactionsFile):
    cgr=CGRUitls()
    [reactions,cntValid,cntAll]=cgr.cgr2SMIRKS(cgrFile)
    cgr.dumpListToFile(reactions,reactionsFile)
        
        
    


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str,
                    help="Input CGR Filename")
    parser.add_argument("outfile",
                    help="Output Reactions  SMILES")
    args = parser.parse_args()
    
    main(args.infile, args.outfile)
    