# %%
import pandas as pd
import numpy as np
import xarray as xa
import matplotlib.pyplot as plt
from combat_atlas import data, config
from combat_atlas.common.caching import compose, lazy, XArrayCache


# %%
class _analysis:
    storage = config.cache/'playground1'
    
    @compose(property, lazy, XArrayCache())
    def pseudobulk(self):
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
            Annotation_minor_subset = [x for (_, x, _), _, _ in x1],
            num_cells = [x for _, x, _ in x1],
        ))
        x3['sample_id'] = x3.scRNASeq_sample_ID + ':' + x3.Annotation_major_subset + ':' + x3.Annotation_minor_subset
        x3 = x3.set_index('sample_id').to_xarray()
        x4 = x.var.reset_index().rename(columns={'index': 'feature_id'})
        x4 = x4.set_index('feature_id').to_xarray()
        x2 = np.stack([x for _, _, x in x1], axis=0)
        x2 = xa.DataArray(x2, [x3.sample_id, x4.feature_id], name='mat')
        x1 = xa.merge([x2, x3, x4])

        x5 = o.columns
        x5 = x5[~x5.str.contains('TCR|BCR|QC|Annotation|Pool_ID|Channel_ID|GEX_region|row')]
        x5 = o[x5.to_list()]

        x6 = x5.drop(columns=['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID', 'Hospitalstay', 'TimeSinceOnset', 'Source'])
        x6 = x6.drop_duplicates().set_index('COMBAT_ID').to_xarray()

        x7 = x5[['COMBAT_participant_timepoint_ID', 'scRNASeq_sample_ID', 'COMBAT_ID', 'Hospitalstay', 'TimeSinceOnset', 'Source']]
        x7 = x7.drop_duplicates().set_index('scRNASeq_sample_ID').to_xarray()

        x8 = xa.merge([
            x1.rename(scRNASeq_sample_ID='sample_scRNASeq_sample_ID'),
            x6,
            x7.rename(COMBAT_ID='sample_COMBAT_ID')
        ])

        return x8

analysis = _analysis()

# %%
if __name__ == '__main__':
    self = analysis

    # %%