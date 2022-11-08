# %%
if __name__ == '__main__':
    __package__ = 'combat_atlas'

# %%
import pandas as pd
import numpy as np
import xarray as xa
from plotnine import *
import statsmodels.api as sm
from .common import helpers
from . import data, config
from .common.caching import compose, lazy, XArrayCache
from .playground1 import analysis as playground1

# %%
x1 = playground1.annot
x2 = x1.cell_scRNASeq_sample_ID
x2 = x1.cells_COMBAT_participant_timepoint_ID.loc[x2]
x2 = x1[['clin_Source', 'clin_COMBAT_ID']].sel(COMBAT_participant_timepoint_ID=x2)
x2 = x2.drop(['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID'])
x1 = xa.merge([x1.cell_Annotation_major_subset, x2])
x1 = x1.sel(cell_id=x1.cell_Annotation_major_subset!='nan')
x1 = x1.to_dataframe()

x2 = list(x1.columns)
x1 = sm.stats.Table.from_data(x1)
x1 = pd.concat([
    v.stack().rename(k)
    for k, v in [
        ('table', x1.table_orig), 
        ('resid', x1.resid_pearson), 
        ('fit', x1.fittedvalues)
    ]
], axis=1).reset_index()
x1['delta'] = x1['table'].astype(str) + '\n' + x1['fit'].astype(int).astype(str)

print(
    ggplot(x1)+
        aes(x2[1], x2[0])+
        geom_tile(aes(fill='resid'))+
        geom_text(aes(label='delta'), size=5)+
        scale_fill_gradient2(
            low='blue', mid='white', high='red',
            midpoint=0
        )+
        labs(x='', y='')+
        theme(
            axis_text_x=element_text(angle=45, hjust=1)
        )
)

# %%
x1 = playground1.annot
x2 = x1.cell_scRNASeq_sample_ID
x2 = x1.cells_COMBAT_participant_timepoint_ID.loc[x2]
x2 = x1[['clin_Source', 'clin_COMBAT_ID']].sel(COMBAT_participant_timepoint_ID=x2)
x2 = x2.drop(['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID'])
x1 = xa.merge([x1.cell_Annotation_major_subset, x2])
x1 = x1.sel(cell_id=x1.cell_Annotation_major_subset!='nan')
x1 = x1.to_dataframe()

x2 = x1.groupby(['cell_Annotation_major_subset', 'clin_Source', 'clin_COMBAT_ID']).size()
x3 = x1.groupby(['clin_Source', 'clin_COMBAT_ID']).size()
x4 = pd.merge(
    x2.rename('n').reset_index(), 
    x3.rename('m').reset_index()
)
x4['freq'] = x4.n/x4.m
x4 = x4[x4['cell_Annotation_major_subset']=='NK']
print(
    ggplot(x4)+aes('clin_Source', '100*freq')+
        geom_violin(aes(fill='clin_Source'))+
        geom_boxplot(aes(fill='clin_Source'), width=0.1)+
        geom_jitter(width=0.1)+
        labs(x='', y='NK cells frequency (%)')+
        theme(
            figure_size=(5, 4),
            legend_position='none',
            axis_text_x=element_text(angle=45, hjust=1)
        )
)
# %%
