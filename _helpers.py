#%%
if __name__ == '__main__':
    __package__ = 'ards_geo'

from pathlib import Path
import xarray as xa
import numpy as np
import pandas as pd
import itertools as it

#%%
class _config:
    project='combat-atlas'
    cache = Path.home()/'.cache'/project
    root = Path(__file__).parent

config = _config()

#%%
def groupby(x1, group_id='group_id', x3 = None, f = None):
    if group_id is None:
        i = 0
        while 'dim_'+str(i) in x1:
            i = i + 1
        group_id = 'dim_'+str(i)
    
    _group_id = '_'+group_id
    x1 = x1.to_dataframe()
    x2 = x1.drop_duplicates().copy()
    x2[group_id] = range(x2.shape[0])
    x1 = x1.reset_index().\
        merge(x2)[[group_id]+x1.index.names].\
        set_index(x1.index.names).to_xarray().\
        rename({group_id: _group_id})
    x2 = x2.set_index(group_id).to_xarray()
    x1 = xa.merge([x1, x2], join='inner')
    if f is None:
        return x1

    x3 = x3.groupby(x1[_group_id]).apply(f)
    x3 = xa.merge([
        x1,
        x3.rename({_group_id: group_id})
    ], join='inner')
    return x3
