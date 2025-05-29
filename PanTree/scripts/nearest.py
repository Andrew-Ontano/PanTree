#!/usr/bin/env python3

import os
import logging
import pandas as pd
import pantreelib as pl
import re
import math
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_vcf):
        logger.error(f"Couldn't find input VCF '{args.input_vcf}. Exiting.'")
        exit()

    # Load VCF
    vcfReport = pl.readVCF(args.input_vcf)
    if args.reference:
        sampleNames = list(vcfReport[0].columns[8:])
        sampleNames[0] = 'Reference'
    else:
        sampleNames = list(vcfReport[0].columns[9:])

    filteredVCF = pl.filterByMissing(vcfReport[0], args.missing)

    targetList = [a.strip() for a in args.targets.split(",")]
#    for record in filteredVCF:
#       pass
    if len([a for a in targetList if a in sampleNames]) != len(targetList):
        logger.error(f"Target names ({','.join(targetList)}) not found in samples ({','.join(sampleNames)}) . Exiting.")
        exit()

    if "Reference" in targetList:
        #alleles = {a: filteredVCF[a] for a in [b for b in targetList if b != "Reference"]}
        #alleles["Reference"] = pd.Series(["0"]*filteredVCF.shape[0])
        #print(alleles)
        targetLength = len(targetList)
        for rowID, rowReport in filteredVCF.iterrows():
            alleles = {a: rowReport[a] for a in [b for b in targetList if b != "Reference"]}
            alleles["Reference"] = 0
            if len(set([alleles[a] for a in alleles])) == targetLength:
                print(rowReport)

    else:
        pass