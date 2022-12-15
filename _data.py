#%%
if __name__ == '__main__':
    __package__ = 'combat_atlas'

import tarfile
import pandas as pd
import xarray as xa
from .common.caching import compose, lazy
from ._helpers import config

#%%
class _data:
    @compose(property, lazy)
    def clinvar(self):
        x = config.cache/'download'/'CBD-KEY-CLINVAR'/'COMBAT_CLINVAR_for_processed.txt'
        x = pd.read_csv(x, sep='\t')
        x = x.set_index('row_number').to_xarray()
        return x

    @compose(property, lazy)
    def rnaseq_wb(self):
        x = config.cache/'download'/'CBD-KEY-RNASEQ-WB'/'Raw_count_data_143_60683.txt'
        x = pd.read_csv(x, sep='\t')
        x = xa.DataArray(
            x, coords=[('feature', x.index), ('sample', x.columns)],
            name='rnaseq'
        )
        return x 

    @compose(property, lazy)
    def rnaseq_wb_logcpm(self):
        x = config.cache/'download'/'CBD-KEY-RNASEQ-WB'/'Logcpm_143_23063.txt'
        x = pd.read_csv(x, sep='\t')
        x = xa.DataArray(
            x, coords=[('feature', x.index), ('sample', x.columns)],
            name='rnaseq'
        )
        return x 

    @compose(property, lazy)
    def citeseq1(self):
        import scanpy
        x = config.cache/'download'/'COMBAT-CITESeq-DATA.h5ad'
        x = scanpy.read(x)
        return x


    def make_cache():
        pass


data = _data()

#%%
if __name__ == '__main__':
    self = data

# %%
    z = self.clinvar
    z = z.sel(row_number=~z.RNASeq_sample_ID.isnull())
    z = z.swap_dims(row_number='RNASeq_sample_ID')
    z = z.rename(RNASeq_sample_ID='sample')

    x2 = self.rnaseq_wb
    x2 = xa.merge([x2, z], join='inner')

# %%
