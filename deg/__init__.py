# %%
if __name__ == '__main__':
    import os
    import sys
    sys.path[0] = os.path.expanduser('~/projects/')
    __package__ = 'combat_atlas.deg'

# %%
import numpy as np
import xarray as xa
import pandas as pd

from ..common import helpers
from ..common.caching import compose, lazy, XArrayCache
from .._helpers import config, groupby

# %%
class _analysis:    
    storage = config.cache/'deg'

    @property
    def data1(self):
        from ..playground1 import analysis as playground1
        x = playground1.pseudobulk1
        x = x.rename(pseudo1_id='group_id1')
        x = x.rename({k: k[8:] for k in x.keys()})
        x = x.rename(
            Annotation_major_subset = 'major',
            Annotation_minor_subset = 'minor',
            mat = 'counts',
            num_cells = 'n'
        )
        return x

    @property
    def data2(self):
        from ..playground1 import analysis as playground1
        x1 = self.data1.copy()
        x2 = playground1.annot.cells_COMBAT_participant_timepoint_ID.sel(
            scRNASeq_sample_ID=x1.scRNASeq_sample_ID
        ).drop('scRNASeq_sample_ID')
        x2 = playground1.annot[['clin_Source', 'clin_COMBAT_ID']].sel(
            COMBAT_participant_timepoint_ID=x2
        ).drop('COMBAT_participant_timepoint_ID')
        x2 = x2.rename(
            clin_Source='group',
            clin_COMBAT_ID='donor'
        )
        x1 = xa.merge([x1, x2])

        x1 = groupby(
            x1[['major', 'group', 'donor']], 'group_id2',
            x1[['counts', 'n']],
            lambda x: x.sum(dim='group_id1')
        )

        x2 = groupby(
            x1[['group', 'donor']], 'group_id3', 
            x1.n, lambda x: x.sum(dim='group_id2')
        )
        x2 = x2.n[x2._group_id3].drop('group_id3')
        x1['freq'] = x1.n/x2
        x1 = x1.drop('n')

        return x1

    @compose(property, lazy, XArrayCache())
    def voom1(self):
        from rpy2.robjects import r as R
        from ..rpy import r_func
        R.setwd(str(config.root))
        R.source('.Rprofile')
        R.source('deg/limma.R')

        x1 = self.data2.copy()
        x1 = x1.drop_dims('group_id1')
        x1['counts'] = x1.counts.transpose('feature_id', 'group_id2')
        x2 = groupby(x1[['major']], 'group_id3')
        x1 = xa.merge([x1[['counts', 'group']], x2], join='inner')
        x6 = x1.group.to_series().drop_duplicates().to_list()
        x6 = pd.CategoricalDtype(x6)
        
        x2 = list()
        #x1 = x1.sel(group_id2=x1._group_id3.isin([0,1,2]))
        for group_id3, x3 in x1.groupby('_group_id3'):
            # group_id3, x3 = next(iter(x1.groupby('_group_id3')))
            print(group_id3)

            x4 = x3.group.to_dataframe()
            x4['group'] = x4.group.astype(x6)

            x5 = r_func(R.fit3)(x3.counts, x4, '~0+group')
            x5 = xa.merge([v.rename(k) for k, v in x5.items()])
            x5['var'] = 'var', x6.categories
            x5['var1'] = 'var1', x5['var'].data
            x5 = x5.expand_dims(group_id3=[group_id3])
            x2.append(x5)
        x2 = xa.concat(x2, 'group_id3')
        x2 = xa.merge([
            x2, 
            x1[['major', '_group_id3']]
        ], join='inner')
        return x2

    @compose(property, lazy, XArrayCache())
    def model2(self):
        from rpy2.robjects import r as R
        from ..rpy import r_func
        R.setwd(str(config.root))
        R.source('.Rprofile')
        R.source('deg/limma.R')

        x1 = self.voom1.copy()
        x2 = x1['var'].to_series().drop_duplicates().to_list()
        x2 = [['HV', x] for x in set(x2)-set(['HV'])]

        i = 0
        x3 = list()
        #x1 = x1.sel(group_id3=[0,1,2])
        for group_id3, x4 in x1.groupby('group_id3'):
            #group_id3, x4 = next(iter(x1.groupby('group_id3')))
            print(group_id3)
            for x5 in x2:
                #x5 = x2[0]
                print(x5)
                x6 = xa.DataArray(
                    np.zeros((x1.sizes['var'],2)),
                    (x1['var'], ('var2', ['control', 'case']))
                )
                x6.loc[x5, :] = [[1,-1],[0,1]]
                x6 = r_func(R.fit4)(
                    *list(x4[['coef', 'cov', 'std', 'sigma', 'df']].values()),
                    x6
                )
                x6 = xa.merge([v.rename(k) for k, v in x6.items()])
                x6 = x6.rename(var2='var')
                x6 = x6.expand_dims(group_id4=[i])
                x6['control'] = 'group_id4', [x5[0]]
                x6['case'] = 'group_id4', [x5[1]]
                x6['_group_id3'] = 'group_id4', [group_id3]
                x3.append(x6)
                i = i + 1
        x3 = xa.concat(x3, 'group_id4')
        x3 = xa.merge([x3, x1[['major']]], join='inner')
        
        return x3

    @compose(property, lazy, XArrayCache())
    def voom2(self):
        from rpy2.robjects import r as R
        from ..rpy import r_func
        R.setwd(str(config.root))
        R.source('.Rprofile')
        R.source('deg/limma.R')

        x1 = self.data2.copy()
        x1 = x1.drop_dims('group_id1')
        x1['counts'] = x1.counts.transpose('feature_id', 'group_id2')

        x2 = groupby(x1[['major']], 'group_id3')
        x1 = xa.merge([x1[['counts', 'group']], x2], join='inner')

        x2 = x1[['_group_id3', 'group']].to_dataframe()
        x2['var'] = x1._group_id3.astype(str) + x2.group
        x2['var'] = x2['var'].astype('category').cat.codes.astype(str)
        x3 = r_func(R.fit3)(x1.counts, x2[['var']], '~0+var')
        x3 = xa.merge([v.rename(k) for k, v in x3.items()])
        x3['var'] = 'var', x3['var'].to_series().str.replace('var', '').astype(int)
        x3['var1'] = 'var1', x3['var1'].to_series().str.replace('var', '').astype(int)
        x2['var'] = x2['var'].astype(int)    
        x2 = x2.drop_duplicates().set_index('var').\
            to_xarray().rename(_group_id3='var_group_id3')
        x4 = xa.merge([
            x3, x2, 
            x1[['_group_id3', 'major']]
        ])
        x4 = x4.sel(var1=x4['var'].data)
        return x4

    @compose(property, lazy, XArrayCache())
    def model3(self):
        from ..rpy import r_func
        from rpy2.robjects import r as R
        R.setwd(str(config.root))
        R.source('.Rprofile')
        R.source('deg/limma.R')

        x1 = self.voom2.copy()
        x2 = x1['group'].to_series().drop_duplicates().to_list()
        x2 = [['HV', x] for x in set(x2)-set(['HV'])]

        i = 0
        x4 = []            
        x3 = x1[['var_group_id3', 'group']].to_dataframe()
        #x3 = x3[x3.var_group_id3.isin([0,1,2])]
        for group_id3, x5 in x3.groupby('var_group_id3'):
            #group_id3, x5 = next(iter(x3.groupby('var_group_id3')))
            print(group_id3)            
            x5 = x5.reset_index().set_index('group')['var']
            for x6 in x2:                
                #x6 = x2[0]
                print(x6)
                x8 = [x5.get(x, None) for x in x6]
                if any(x is None for x in x8):
                    continue

                x7 = xa.DataArray(
                    np.zeros((x1.sizes['var'], 2)),
                    [x1['var'], ('var2', ['control', 'case'])]
                )
                x7.loc[x8, :] = [[1,-1], [0,1]]
                x7 = r_func(R.fit4)(
                    *list(x1[['coef', 'cov', 'std', 'sigma', 'df']].values()),
                    x7
                )
                x7 = xa.merge([v.rename(k) for k, v in x7.items()])
                x7 = x7.rename(var2='var')
                x7 = x7.expand_dims(group_id4=[i])
                x7['control'] = 'group_id4', [x6[0]]
                x7['case'] = 'group_id4', [x6[1]]
                x7['_group_id3'] = 'group_id4', [group_id3]
                x4.append(x7)
                i = i + 1        
        x4 = xa.concat(x4, dim='group_id4')
        x4 = xa.merge([x4, x1[['major']]], join='inner')
        
        return x4

    @compose(property, lazy, XArrayCache())
    def samples2(self):
        x = self.data2.drop_dims(['feature_id', 'group_id1'])
        return x

    @compose(property, lazy, XArrayCache())
    def genes2(self):
        return self.data2.feature_id.rename('symbol')

    def rnaseq_counts2(self, purification, cell_type):  
        analysis = self
        class counts:      
            storage = analysis.storage/'rnaseq_counts2'/(purification+'_'+cell_type)

            @compose(property, lazy, XArrayCache())
            def data(self):
                import sparse
                x = analysis.data2.drop_dims('group_id1')
                x = x.sel(group_id2=x['purification']==purification)
                x = x.sel(
                    group_id2=x['blueprint.labels']==cell_type
                )
                x = x['counts']
                x.data = sparse.COO(x.data)
                return x

        counts = counts()

        return counts.data.todense()

analysis = _analysis()

# %%
if __name__ == '__main__':
    self = analysis   

    # %%