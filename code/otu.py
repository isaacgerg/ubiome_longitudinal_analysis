import json
import datetime
from collections import defaultdict

import numpy as np

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
        
        # form Tree
        for k in lst:
            self._tree[str(k['taxon'])] = (k['tax_rank'], k['tax_name'], k['count_norm']/1000000, str(k['parent']))        
        
    def getTaxonomy(self, taxon='genus'):   
        taxNames = []        
        taxCounts = defaultdict(lambda: 0)
        for k in self._tree.keys():            
            if self._tree[k][0] == taxon:
                #taxNames.append(self._tree[k][2])
                taxCounts[self._tree[k][1]] += self._tree[k][2]

        taxNames = list(taxCounts.keys())
    
        N = len(taxNames)
        r = np.zeros(N)
        for taxName,idx in zip(taxNames,range(N)):
            r[idx] = taxCounts[taxName]
            
        return (taxNames, r)