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
from .common.caching import compose, lazy, XArrayCache, CSVCache
from .playground1 import analysis as playground1

storage = config.cache/'playground2'

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
from types import SimpleNamespace

class _analysis1:
    @compose(property, lazy, XArrayCache())
    def gsva(self):
        from .sigs._sigs import sigs
        from .sigs.gsva import gsva
        x2 = self.x2
        x1 = sigs.all1
        x1 = x1.sel(sig=x1.sig_prefix.isin(["KEGG1", "HALLMARK"]))
        x1 = xa.merge([
            x1,
            self.feature_entrez.rename(entrez='gene')
        ], join='inner')
        x1 = xa.merge([
            x1,
            self.log1p_rpm.log1p_rpm.rename(group='sample')
        ], join='inner')
        x1['log1p_rpm'] = xa.apply_ufunc(
            np.matmul, x1.feature_entrez, x1.log1p_rpm,
            input_core_dims=[['gene', x2.feature], [x2.feature, 'sample']],
            output_core_dims=[['gene', 'sample']]
        )
        x1['log1p_rpm'] = x1.log1p_rpm/x1.feature_entrez.sum(dim=x2.feature)
        x1['gene'] = x1.gene.astype(np.int32)

        x3 = gsva(
            x1.log1p_rpm,
            x1.set.to_series_sparse().reset_index().drop(columns=['set'])
        )
        x3 = xa.merge([x3.rename('gsva').to_dataset(), x1[['sig_prefix']]])
        x3['sample'] = x3['sample'].astype(x1.sample.dtype)
        x3 = x3.rename(sample='group')
        return x3.gsva
        
    @compose(property, lazy)
    def summary1(self):
        x2 = self.x2
        s = xa.merge([
            self.gsva,
            self.log1p_rpm[x2.diag]
        ])
        s = s.groupby(x2.diag).apply(lambda x: xa.Dataset(dict(
            mu = x.gsva.mean(dim='group'),
            sigma = x.gsva.std(dim='group'),
            n = x.sizes['group']
        )))
        return s

    @compose(property, lazy, CSVCache(ext='.csv'))
    def summary2(self):
        x7 = self.x2.summary2
        x2 = self.summary1
        x3 = x2.sel({self.x2.diag: x7[0]})
        x4 = x2.sel({self.x2.diag: x7[1]})

        x5 = x4.mu-x3.mu
        x6 = (x4.n*(x4.sigma**2)+x3.n*(x3.sigma**2))/(x4.n+x3.n)
        x5 = x5/x6
        x5 = xa.merge([x5.to_dataset(self.x2.diag), (x3.mu/x3.sigma).rename(x7[0])])

        x5 = x5.to_dataframe().reset_index()
        x5 = x5[['sig','sig_prefix',x7[0]]+x7[1]]
        x5['score'] = x5[x7[1]].abs().max(axis=1)
        x5['sig'] = x5.sig.str.replace('^[^_]*_', '', regex=True)
        x5 = x5.sort_values('score', ascending=False)
        return x5

class _combat_analysis1(_analysis1):
    storage = storage/'analysis1'

    x2 = SimpleNamespace(
        diag='clin_Source',
        cells='cell_type',
        NK='NK',
        counts='pseudo_mat',
        donor='clin_COMBAT_ID',
        cell='group',
        feature='feature_id',
        summary2=['HV', ['COVID_CRIT', 'COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'Sepsis']]
    )

    @property
    def feature_entrez(self):
        return playground1.feature_entrez.rename(entrez_id='entrez').rename('feature_entrez')

    @property
    def data1(self):
        x2 = self.x2
        x1 = playground1.data1
        x1 = x1.swap_dims(gene_id='feature_id')
        x1 = x1[[x2.counts, x2.donor, x2.diag, x2.cells]]
        x1 = x1.stack({'group1': ['sample_id', x2.cells]})
        x1 = x1.assign_coords(group=('group1', x1['group1'].to_series().astype('category').cat.codes))
        x1['sample_id1'] = 'group', x1.sample_id.data
        x1['cell_type1'] = 'group', x1.cell_type.data
        x1 = x1.swap_dims(group1='group').drop('group1')
        x1 = x1.rename(sample_id1='sample_id', cell_type1=x2.cells)
        x1 = x1.sel({'group': x1[x2.cells]==x2.NK})
        x1 = x1.drop('gene_id')
        return x1

    @property
    def log1p_rpm(self):
        x1 = self.data1
        x2 = self.x2

        x1['group1'] = (x1[x2.donor] + x1[x2.cells]).to_series().\
            astype('category').cat.codes.to_xarray()
        x3 = x1.groupby('group1')
        x3 = x3.apply(lambda x: xa.merge([
            x[[x2.donor, x2.diag, x2.cells]].isel(group=0),
            x[x2.counts].mean(dim='group')
        ]))
        x3['total'] = x3[x2.counts].sum(dim=x2.feature)
        x3['log1p_rpm'] = np.log1p(x3[x2.counts]/x3['total'])
        x3 = x3.drop('group').rename(group1='group')
        x3 = x3.sel(group=x3.total>0)
        return x3

x1 = _combat_analysis1()

self = x1
# %%
