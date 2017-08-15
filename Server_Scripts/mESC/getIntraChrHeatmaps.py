from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome
from mirnylib import h5dict

def extractResolutionFromFileName(fname):
    try:
        raw_heatmap = h5dict.h5dict(fname, mode='r') #open heatmap
        resolution = int(raw_heatmap['resolution']) #get the resolution
        del raw_heatmap #close heatmap
        return resolution
    except:
        try:
            if "/" in fname:
                fname = fname.split("/")[-1]
            res = fname.split("res")[-1].split("k")[0]
            res = int(res)*1000
            return res
        except:
            print "Warning! Unable to resolve resolution from file name"
            return None


def get_chromosomes(hm_file,genome_db,resolution,chrNumb=None):
        if extractResolutionFromFileName(hm_file) != resolution:
                print "WARNING! Provided resolution ",resolution, "does not match ", extractResolutionFromFileName(hm_file),"extracted from file name ",hm_file
        if "hiRes.hm" in hm_file:
                type = "HiRes"
        elif "bychr.hm" in hm_file:
                type = "bychr"
        else:
            print "Warning: cannot resolve type of data from filename"
            try:
                print "Warning: trying hires hic"
                raw_heatmap = h5dict.h5dict(fname, mode='r') #open heatmap
                if "0 0" in raw_heatmap.keys():
                    type = "HiRes"
                else:
                    print "HiRes hic Failed! Assuming bychr type"
                    type = "bychr"
            except:
                print "HiRes hic Failed! Assuming bychr type"
                type = "bychr"
        if type=="HiRes":
                from hiclib import highResBinnedData
                # Create a  object, load the data.
                print "creating an object"
                hmap = highResBinnedData.HiResHiC(genome_db,resolution)
                print "loading data"
                hmap.loadData(hm_file, mode="cis")
                print "Data loaded"
                if chrNumb != None:
                    return hmap.data[(chrNumb,chrNumb)].getData()
                return [hmap.data[(i,i)].getData() for i in xrange(genome_db.chrmCount)]
                #cisKeys are tuples like (N,N) where N is 0..Number_of_chrms-1
        elif type=="bychr":
                from hiclib import binnedData
                print "creating an object"
                hmap = binnedData.binnedData(resolution,genome_db)

                print "loading data"
                hmap.simpleLoad(hm_file,"heatmap")
                data=hmap.dataDict["heatmap"]
                assert len(data)==genome_db.numBins
                print "Data loaded"
                if chrNumb != None:
                	return data[genome_db.chrmStartsBinCont[chrNumb]:genome_db.chrmEndsBinCont[chrNumb],genome_db.chrmStartsBinCont[chrNumb]:genome_db.chrmEndsBinCont[chrNumb]]
                return [data[genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i],genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i]] for i in xrange(genome_db.chrmCount)]
        else:
                raise "Error: can not recognize heatmap format from file name"