# %%
import pandas as pd
import numpy as np
import xarray as xa
import matplotlib.pyplot as plt
from plotnine import *
import combat_atlas.common.helpers
from combat_atlas import data, config
from combat_atlas.common.caching import compose, lazy, XArrayCache


# %%
class _analysis:
    storage = config.cache/'playground1'

    @compose(property, lazy, XArrayCache())
    def annot(self):
        x = data.citeseq1

        o = x.obs.copy()

        x1 = o[['scRNASeq_sample_ID', 'Annotation_major_subset', 'Annotation_minor_subset']]
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

        o = x.obs.copy()   
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

analysis = _analysis()

# %%
if __name__ == '__main__':
    self = analysis

    # %%
    x1 = self.annot.feature_gene_ids.sel(feature_id='C5AR1')

    x2 = self.bulk
    x2 = 1e6*x2/x2.sum(dim='gene_id')
    x2 = x2.sel(gene_id=x1.data)
    x2 = xa.merge([
        x2,
        self.annot[['RNASeq_sample_ID', 'bulk_COMBAT_participant_timepoint_ID']].\
                rename(bulk_COMBAT_participant_timepoint_ID='COMBAT_participant_timepoint_ID')
    ]).to_dataframe().reset_index()

    x3 = self.pseudobulk2.pseudo2_mat
    x3 = 1e6*x3/x3.sum(dim='feature_id')
    x3 = x3.sel(feature_id=x1.feature_id.data)
    x3 = xa.merge([
        x3.rename(pseudo2_id='scRNASeq_sample_ID'),
        self.annot[['scRNASeq_sample_ID', 'cells_COMBAT_participant_timepoint_ID']].\
                rename(cells_COMBAT_participant_timepoint_ID='COMBAT_participant_timepoint_ID')                
    ]).to_dataframe().reset_index()

    x4 = pd.merge(x2, x3)

    print(
        ggplot(x4)+aes('bulk_mat', 'pseudo2_mat')+
            geom_point()+
            geom_smooth(method='lm')+
            geom_hline(yintercept=x4.pseudo2_mat.mean(), linetype='--')+
            labs(
                x='bulk (RPM)', y='pseudbulkd (RPM)', 
                title=f'C5AR1 (R={x4[["bulk_mat", "pseudo2_mat"]].corr().to_numpy()[0,1]:.2f})'
            )
    )

    # %%
    x1 = self.annot.feature_gene_ids.sel(feature_id='C5AR1')

    x2 = self.bulk
    x2 = 1e6*x2/x2.sum(dim='gene_id')
    x2 = x2.sel(gene_id=x1)
    x2 = xa.merge([
        x2,
        self.annot[['RNASeq_sample_ID', 'bulk_COMBAT_participant_timepoint_ID']].\
            rename(bulk_COMBAT_participant_timepoint_ID='COMBAT_participant_timepoint_ID')
    ]).to_dataframe()
    x2 = x2.merge(
        self.annot.drop_dims(list(set(self.annot.dims)-set(['COMBAT_participant_timepoint_ID']))).\
            rename(clin_COMBAT_ID='COMBAT_ID').to_dataframe().reset_index()
    )
    x2 = x2.merge(
        self.annot.drop_dims(list(set(self.annot.dims)-set(['COMBAT_ID']))).\
            to_dataframe().reset_index()
    )

    print(
        ggplot(x2)+aes('clin_Source', 'bulk_mat')+
            geom_violin(aes(fill='clin_Source'))+
            geom_jitter(width=0.2)+
            coord_flip()+
            labs(x='C5AR1 bulk (RPM)')+
            theme(legend_position='none')
    )

    # %%
    x1 = self.annot.feature_gene_ids.sel(feature_id='C5AR1')

    x2 = self.pseudobulk2
    x2 = 1e6*x2/x2.sum(dim='feature_id')
    x2 = x2.sel(feature_id=x1.feature_id)
    x2 = xa.merge([
        x2,
        self.annot[['scRNASeq_sample_ID', 'cells_COMBAT_participant_timepoint_ID']].\
            rename(
                scRNASeq_sample_ID='pseudo2_id',
                cells_COMBAT_participant_timepoint_ID='COMBAT_participant_timepoint_ID'
            )
    ]).to_dataframe()
    x2 = x2.merge(
        self.annot.drop_dims(list(set(self.annot.dims)-set(['COMBAT_participant_timepoint_ID']))).\
            rename(clin_COMBAT_ID='COMBAT_ID').to_dataframe().reset_index()
    )
    x2 = x2.merge(
        self.annot.drop_dims(list(set(self.annot.dims)-set(['COMBAT_ID']))).\
            to_dataframe().reset_index()
    )

    print(
        ggplot(x2)+aes('clin_Source', 'pseudo2_mat')+
            geom_violin(aes(fill='clin_Source'))+
            geom_jitter(width=0.2)+
            coord_flip()+
            labs(x='C5AR1 pseudobulk (RPM)')+
            theme(legend_position='none')
    )

    
# %%
