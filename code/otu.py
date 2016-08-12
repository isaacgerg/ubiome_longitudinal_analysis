import json
import datetime

def tree(): return defaultdict(tree)

class OTU:
    def __init__(self, jsonFilename):
        self._filename = jsonFilename
        r = open(self._filename).read()
        self._json = json.loads(r)
        k = self._json['sampling_time']
        self._sampleTime = datetime.datetime(int(k[0:4]), int(k[5:7]), int(k[8:10]))

        self._description = 'TODO'
        
        self._tree = tree()
        
        lst = self._json['ubiome_bacteriacounts']
        
        self._phyla = []
        
    def getPhyla(self):
        if self._phyla is not []:
            return self._phyla
           
        lst = self._json['ubiome_bacteriacounts']
        for k in lst:
            if k['tax_rank'] == 'phylum':
                self._phyla.append(k['tax_name'])

        self._phyla = list(set(self._phyla))
        
        #return self._phyla
    
        # Get counts for each phyla from each dataset
        N = len(fn)
        p = {}
        for k in phyla:
            p[k] = np.zeros(N)
            c = 0
            for f in fn:
                r = open(f).read()
                lst = json.loads(r)['ubiome_bacteriacounts'] 
                for kk in lst:
                    if kk['tax_rank'] == 'phylum' and kk['tax_name'] == k:
                        p[kk['tax_name']][c] =  kk['count_norm']/10000
                c += 1    