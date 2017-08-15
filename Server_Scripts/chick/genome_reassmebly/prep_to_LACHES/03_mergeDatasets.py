from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome 
import numpy as np 
import os
import sys

genomeName = "GalGal5filtered"
genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")

basefolder = "/mnt/storage/home/vsfishman/HiC/data/chick/mapped-GalGal5filtered/B1_TTAGGC_L001_/"
filename = "chunk0001.hdf5"
				
TR = HiCdataset(basefolder+filename+".HiCdataset", genome=genome_db,
                                    maximumMoleculeLength=500,enzymeName = "HindIII",tmpFolder = "tmp",
                                    mode='w')  # remove inMemory if you don't have enough RAM
TR.parseInputData(dictLike=basefolder+filename)
TR.filterDuplicates()
TR.filterLarge(10000,10)
TR.filterExtreme(cutH=0.001, cutL=0)
TR.writeFilteringStats()
TR.printMetadata(saveTo=basefolder+filename+".stat")
TR.saveHeatmap(basefolder+filename+".hm-res-1000kb",1000000)
comment ="""
        #------------------------End set of filters applied----------

    print("----->Building Raw heatmap at different resolutions")
    TR.printStats()
    for res in wholeGenomeResolutionsKb:    
        TR.saveHeatmap(out_file + "-{0}k.hm".format(res), res*1000)
    for res in byChromosomeResolutionsKb: 
        TR.saveByChromosomeHeatmap(out_file + "-{0}k.byChr".format(res), res*1000)
    for res in HiResWithOverlapResolutionsKb[:-skip]:
        TR.saveHiResHeatmapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res*1000)        
    for res in SuperHiResWithOverlapResolutionsKb[:-skip]:
        TR.saveSuperHighResMapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res*1000)            



#Now merging different experiments alltogether
#note that the first column is not here, as it is a replica 
experiments = set([(i[0], i[2], i[3]) for i in combinedExperimentNames])
print(experiments)

for experiment in experiments:
    workingGenome = experiment[1]
    myExperimentNames = [i[1] + "_refined.frag" for i in combinedExperimentNames if (i[0], i[2], i[3]) == (experiment[0], experiment[1],experiment[2])]    
    assert len(myExperimentNames) > 0
    if len(myExperimentNames) > 0:
        #If we have more than one experiment (replica) for the same data, we can combine. 
        if genomeName != workingGenome: 
            raise Exception ("Check Genome!")
        TR = HiCdataset(os.path.join(workingGenome, "%s-all-%s_refined.frag" %
                                     (experiment[0],experiment[2])), genome=genome_db,
                                     enzymeName = experiment[2],tmpFolder = "tmp",dictToStoreIDs="h5dict")
        statSaveName = os.path.join("statistics", workingGenome, "%s-all-%s_refined.stat" % (experiment[0], experiment[2]))

        TR.merge(myExperimentNames)
        TR.printMetadata(saveTo=statSaveName)
        for res in wholeGenomeResolutionsKb:    
            TR.saveHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.hm" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in byChromosomeResolutionsKb: 
            TR.saveByChromosomeHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in HiResWithOverlapResolutionsKb:
            TR.saveHiResHeatmapWithOverlaps(os.path.join(workingGenome, "%s-all-%s-{0}k_HighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
        for res in SuperHiResWithOverlapResolutionsKb:
            TR.saveSuperHighResMapWithOverlaps(os.path.join(workingGenome, "%s-all-%s-{0}k_SuperHighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)

"""