# %%
if __name__ == '__main__':
    __package__ = 'combat_atlas.playground1'

# %%
import pandas as pd
import numpy as np
import xarray as xa
from plotnine import *
from ..common import helpers
from .. import data, config
from ..common.caching import compose, lazy, XArrayCache

# %%
class _analysis:
    storage = config.cache/'playground1'

    @compose(property, lazy, XArrayCache())
    def annot(self):
        x = data.citeseq1

        o = x.obs.copy()

        x1 = o[['scRNASeq_sample_ID', 'Annotation_major_subset', 'Annotation_minor_subset']].copy()
        x1['Annotation_major_subset'] = np.where(
            x1.Annotation_minor_subset=='cMono.PLT', 
            'cMono', 
            x1.Annotation_major_subset
        )
        x1['Annotation_minor_subset'] = np.where(
            (x1.Annotation_major_subset!='nan') & (x1.Annotation_minor_subset=='nan'),
            x1.Annotation_major_subset,
            x1.Annotation_minor_subset
        )
        x1 = x1.to_xarray()
        x1 = x1.rename(index='cell_id').rename({k: 'cell_'+k for k in x1.keys()})

        x6 = o[['COMBAT_participant_timepoint_ID', 'COMBAT_ID', 'Hospitalstay', 'TimeSinceOnset', 'Source']]
        x6 = x6.drop_duplicates().set_index('COMBAT_participant_timepoint_ID').to_xarray()
        x6 = x6.rename({k: 'clin_'+k for k in x6.keys()})

        x2 = o[['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID']]
        x2 = x2.drop_duplicates().set_index('scRNASeq_sample_ID').to_xarray()
        x2 = x2.rename({k: 'cells_'+k for k in x2.keys()})

        x3 = o.columns
        x3 = x3[~x3.str.contains('TCR|BCR|QC|Pool_ID|Annotation|Channel_ID|GEX_region|row')]
        x3 = o[x3.to_list()]
        x3 = x3.drop(columns=['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID', 'Hospitalstay', 'TimeSinceOnset', 'Source'])
        x3 = x3.drop_duplicates().set_index('COMBAT_ID').to_xarray()
        x3 = x3.rename({k: 'patient_'+k for k in x3.keys()})

        x5 = x.var.reset_index().rename(columns={'index': 'feature_id'})
        x5 = x5.set_index('feature_id').to_xarray()
        x5 = x5.rename(feature_types='types')
        x5 = x5.rename({k: 'feature_'+k for k in x5.keys()})

        x7 = data.clinvar[['RNASeq_sample_ID', 'COMBAT_participant_timepoint_ID']].to_dataframe()
        x7 = x7[~x7.RNASeq_sample_ID.isna()]
        x7 = x7.set_index('RNASeq_sample_ID').to_xarray()
        x7 = x7.rename({k: 'bulk_'+k for k in x7.keys()})

        x4 = xa.merge([x1, x6, x2, x3, x5, x7], join='inner')
        return x4

    @compose(property, lazy)
    def cells(self):
        import sparse
        x = data.citeseq1
        x1 = xa.DataArray(
            sparse.GCXS(x.layers['raw'], compressed_axes=[0]), 
            [self.annot.cell_id, self.annot.feature_id],
            name='cell_mat'
        )
        return x1

    @compose(property, lazy)
    def pseudobulk3(self):
        #this is the same as psuedobulk1, but much slower than pseudobulk1. not used,
        x = self.cells
        x = xa.merge([
            x,
            self.annot.drop_dims(set(self.annot.dims)-set(['cell_id']))
        ])
        x['pseudo1_id'] = x.cell_scRNASeq_sample_ID + ':' + x.cell_Annotation_major_subset + ':' + x.cell_Annotation_minor_subset
        x = x.groupby('pseudo1_id')
        x = x.apply(lambda x: xa.merge([
            x.cell_mat.todense().sum(dim='cell_id'),
            xa.DataArray(x.cell_mat.sizes['cell_id'], name='cell_num_cells'),
            x[[
                'cell_scRNASeq_sample_ID', 
                'cell_Annotation_major_subset', 
                'cell_Annotation_minor_subset'
            ]].isel(cell_id=0).drop('cell_id')
        ]))
        x = x.rename({
            k: k.replace('cell_', 'pseudo1_') for k in x.keys()
        })
        return x
    
    @compose(property, lazy, XArrayCache())
    def pseudobulk1(self):
        x = data.citeseq1
        #o = x.obs.copy()   

        o = self.annot[[
            'cell_scRNASeq_sample_ID', 
            'cell_Annotation_major_subset', 
            'cell_Annotation_minor_subset'
        ]].sel(cell_id=x.obs.index.to_list())
        o = o.rename({k: k[5:] for k in o.keys()})
        o = o.to_dataframe()        
        o['row'] = range(o.shape[0])
        
        x1 = x.layers['raw']
        x1 = [        
            (i, x1[r.to_list(),])
            for i, r in o.groupby([
                'scRNASeq_sample_ID', 
                'Annotation_major_subset',
                'Annotation_minor_subset'
            ]).row
        ]
        x1 = [
            (i, x.shape[0], np.array(x.todense()).sum(axis=0))
            for i, x in x1
        ]
        x3 = pd.DataFrame(dict(
            scRNASeq_sample_ID = [x for (x, _, _), _, _ in x1],
            Annotation_major_subset = [x for (_, x, _), _, _ in x1],
            Annotation_minor_subset = [x for (_, _, x), _, _ in x1],
            num_cells = [x for _, x, _ in x1],
        ))
        x3['pseudo1_id'] = x3.scRNASeq_sample_ID + ':' + x3.Annotation_major_subset + ':' + x3.Annotation_minor_subset
        x3 = x3.set_index('pseudo1_id').to_xarray()
        x2 = np.stack([x for _, _, x in x1], axis=0)
        x2 = xa.DataArray(x2, [x3.pseudo1_id, self.annot.feature_id], name='mat')
        x1 = xa.merge([x2, x3])
        x1 = x1.rename({k: 'pseudo1_'+k for k in x1.keys()})

        return x1

    @compose(property, lazy)
    def pseudobulk2(self):
        x = self.pseudobulk1.copy()
        x['pseudo2_id'] = x.pseudo1_scRNASeq_sample_ID
        x = x.groupby('pseudo2_id')
        x = x.apply(lambda x: xa.merge([
            x.pseudo1_mat.sum(dim='pseudo1_id').rename('pseudo2_mat'),
            x.pseudo1_num_cells.sum(dim='pseudo1_id').rename('pseudo2_num_cells')
        ]))
        return x

    @compose(property, lazy)
    def bulk(self):
        x = data.rnaseq_wb
        x = x.rename(
            feature='gene_id',
            sample='RNASeq_sample_ID'
        ).rename('bulk_mat')
        return x

    @compose(property, lazy)
    def bulk1(self):
        x = data.rnaseq_wb_logcpm
        x = x.rename(
            feature='gene_id',
            sample='RNASeq_sample_ID'
        ).rename('bulk1_mat')
        return x

    @compose(property, lazy)
    def cell_type_freq1(self):
        import statsmodels.api as sm

        x = self.data1.copy()
        x1 = x[['pseudo_mat', 'bulk_mat', 'pseudo_num_cells', 'feature_id']].copy()
        x1['bulk_mat'] = 1e6*x1.bulk_mat/x1.bulk_mat.sum(dim='gene_id')
        x1['pseudo_mat'] = 1e6*x1.pseudo_mat/x1.pseudo_mat.sum(dim='gene_id')
        x2 = x1.bulk_mat
        x1['bulk_mean'] = x2.mean(dim='sample_id')
        x1['bulk_std'] = x2.std(dim='sample_id')
        x1['bulk_mat'] = (x2-x1.bulk_mean)/x1.bulk_std
        x2 = x1.pseudo_mat
        x1['pseudo_mat'] = (x2-x2.mean(dim='sample_id'))/x2.std(dim='sample_id')
        x1 = x1.sel(gene_id=x1.bulk_mat.isnull().sum(dim='sample_id')==0)
        x1 = x1.fillna(0)
        x2 = xa.apply_ufunc(
            lambda X, y: sm.OLS(y, X).fit().params,
            x1.pseudo_mat, x1.bulk_mat,
            input_core_dims=[['gene_id', 'cell_type'], ['gene_id']],
            output_core_dims=[['cell_type']],
            vectorize=True
        )
        x1['coef'] = x2

        return x1

    @compose(property, lazy)
    def data1(self):
        x1 = self.pseudobulk1
        x1 = x1.rename(
            pseudo1_Annotation_major_subset='cell_type'
        )
        x1 = xa.merge([
            x1, 
            self.annot.feature_gene_ids.rename('gene_id')
        ]).set_coords('gene_id').\
            swap_dims(feature_id='gene_id')

        x3 = self.annot.cells_COMBAT_participant_timepoint_ID.sel(
            scRNASeq_sample_ID=x1.pseudo1_scRNASeq_sample_ID
        ).rename('COMBAT_participant_timepoint_ID').drop('scRNASeq_sample_ID')
        x1 = xa.merge([x1, x3])

        x3 = [
            x1.COMBAT_participant_timepoint_ID.to_series(), 
            x1.cell_type.to_series()
        ]
        x3 = pd.MultiIndex.from_tuples(list(zip(*x3)))
        x3 = xa.DataArray(x3, [x1.pseudo1_id], name='new_id')
        x1 = x1.groupby(x3).apply(lambda x: xa.merge([
            x.pseudo1_mat.sum(dim='pseudo1_id').rename('pseudo_mat'),
            x.pseudo1_num_cells.sum(dim='pseudo1_id').rename('pseudo_num_cells')
        ])).unstack('new_id').rename(
            new_id_level_0='COMBAT_participant_timepoint_ID',
            new_id_level_1='cell_type'
        )
        x1 = x1.fillna(0)
        
        x2 = 2**self.bulk1.rename('bulk_mat')
        x2 = xa.merge([
            x2,
            self.annot.bulk_COMBAT_participant_timepoint_ID.\
                rename('COMBAT_participant_timepoint_ID')
        ], join='inner').set_coords('COMBAT_participant_timepoint_ID').\
            swap_dims(RNASeq_sample_ID='COMBAT_participant_timepoint_ID').\
            drop('RNASeq_sample_ID')

        x1 = xa.merge([x1, x2], join='inner')

        x1 = x1.merge(
            self.annot.drop_dims(set(self.annot.dims)-set(['COMBAT_participant_timepoint_ID'])),
            join='inner'
        )
        x1 = x1.merge(
            self.annot.drop_dims(set(self.annot.dims)-set(['COMBAT_ID'])).\
                sel(COMBAT_ID=x1.clin_COMBAT_ID).drop('COMBAT_ID'),
            join='inner'
        )

        x1 = x1.rename(COMBAT_participant_timepoint_ID='sample_id')

        return x1

    @compose(property, lazy, XArrayCache())
    def feature_entrez(self):
        from ..sigs.entrez import symbol_entrez

        x = self.annot.feature_id.to_series()
        x = symbol_entrez(x)
        x = x.rename(
            symbol='feature_id',
            Entrez_Gene_ID = 'entrez_id'
        )
        return x

    @compose(property, lazy)
    def data2(self):
        import sparse
        from scipy.stats import rankdata
        from ..sigs import sigs
        from ..sigs.fit import fit_gsea
        
        x3 = self.data1.copy()
        x3['pseudo_mat'] = x3.pseudo_mat.sum(dim='cell_type')
        x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
        x3 = x3.drop_dims('cell_type')
        #x3 = x3.sel(sample_id=x3.clin_Source.isin(['HV', 'COVID_MILD']))
        x3 = x3.sel(sample_id=x3.clin_Source.isin(['COVID_CRIT', 'COVID_SEV', 'Sepsis']))

        x2 = xa.apply_ufunc(
            lambda x, y: np.corrcoef(rankdata(x), rankdata(y))[0,1],
            x3.pseudo_mat, x3.bulk_mat,
            input_core_dims=[['sample_id'], ['sample_id']],
            output_core_dims=[[]],
            vectorize=True
        ).rename('R').to_dataset().dropna('gene_id')
        x2 = x2.swap_dims(gene_id='feature_id')
        x2 = xa.merge([x2, self.feature_entrez], join='inner')
        x2['R'] = xa.apply_ufunc(
            np.matmul, x2.symbol_entrez, x2.R,
            input_core_dims=[['entrez_id', 'feature_id'], ['feature_id']],
            output_core_dims=[['entrez_id']]
        )
        x2['R'] = x2.R/x2.symbol_entrez.sum(dim='feature_id')
        x2 = x2.sel(entrez_id=x2.R.dropna('entrez_id').entrez_id)
        x2 = x2.rename(entrez_id='gene')
        x4 = x2.symbol_entrez.to_series_sparse().reset_index()
        x4 = x4.groupby('gene').feature_id.apply(lambda x: ','.join(x)).rename('symbol')
        x4 = x4.to_xarray()
        x2['symbol'] = x4
        x2 = x2.drop_dims('feature_id')

        x5 = x2.R.copy().rename('t')
        x5 = np.abs(x5)
        x5.data = sparse.COO(x5.data)
        x5 = xa.merge([x5, sigs.all1.rename('s')], join='inner')
        x5 = fit_gsea(x5.t.expand_dims(g=[1]), x5.s, 1e5)
        x5 = x5.squeeze('g').drop('g')
        x2 = xa.merge([x2, x5], join='inner')

        return x2

    @compose(property, lazy)
    def data3(self):
        import sparse
        from scipy.stats import rankdata
        from ..sigs import sigs
        from ..sigs.fit import fit_gsea
        
        x3 = self.data1.copy()
        x3 = x3.drop_dims('cell_type')
        x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')
        x4 = {
            'HV': 0,
            'COVID_MILD': 1,
            'COVID_SEV': 2,
            'COVID_CRIT': 3,
            'Sepsis': 3,
        }
        x3['severity'] = 'sample_id', [x4[s] for s in x3.clin_Source.data]

        x2 = xa.apply_ufunc(
            lambda x, y: np.corrcoef(rankdata(x), (y))[0,1],
            x3.bulk_mat, x3.severity,
            input_core_dims=[['sample_id'], ['sample_id']],
            output_core_dims=[[]],
            vectorize=True
        ).rename('R').to_dataset().dropna('gene_id')
        x2 = x2.swap_dims(gene_id='feature_id')
        x2 = xa.merge([x2, self.feature_entrez], join='inner')
        x2['R'] = xa.apply_ufunc(
            np.matmul, x2.symbol_entrez, x2.R,
            input_core_dims=[['entrez_id', 'feature_id'], ['feature_id']],
            output_core_dims=[['entrez_id']]
        )
        x2['R'] = x2.R/x2.symbol_entrez.sum(dim='feature_id')
        x2 = x2.sel(entrez_id=x2.R.dropna('entrez_id').entrez_id)
        x2 = x2.rename(entrez_id='gene')
        x4 = x2.symbol_entrez.to_series_sparse().reset_index()
        x4 = x4.groupby('gene').feature_id.apply(lambda x: ','.join(x)).rename('symbol')
        x4 = x4.to_xarray()
        x2['symbol'] = x4
        x2 = x2.drop_dims('feature_id')

        x5 = x2.R.copy().rename('t')
        x5.data = sparse.COO(x5.data)
        x5 = xa.merge([x5, sigs.all1.rename('s')], join='inner')
        x5 = fit_gsea(x5.t.expand_dims(g=[1]), x5.s, 1e5)
        x5 = x5.squeeze('g').drop('g')
        x2 = xa.merge([x2, x5], join='inner')

        return x2



analysis = _analysis()

# %%
if __name__ == '__main__':
    self = analysis

# %%