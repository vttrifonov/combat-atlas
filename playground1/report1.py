# %%
from plotnine import *
import matplotlib.pyplot as plt
import numpy as np
from combat_atlas.playground1 import analysis

# %%
self = analysis

# %%
x3 = self.data1.copy()
x3['bulk_mat'] = 1e6*x3.bulk_mat/x3.bulk_mat.sum(dim='gene_id')
x3['pseudo_mat'] = x3.pseudo_mat.sum(dim='cell_type')
x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
x3 = x3.drop_dims('cell_type')
x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index()    
print(
    ggplot(x2)+aes('bulk_mat', 'pseudo_mat')+
        geom_point()+
        geom_smooth(method='lm')+
        geom_hline(yintercept=x2.pseudo_mat.mean(), linetype='--')+
        labs(
            x='bulk (RPM)', y='pseudbulkd (RPM)', 
            title=f'{x2.feature_id.iloc[0]} (R={x2[["bulk_mat", "pseudo_mat"]].corr().to_numpy()[0,1]:.2f})'
        )
)

# %%
x2 = x3.isel(sample_id=0).to_dataframe().reset_index()
x2['f'] = x2.feature_id.str.contains("^RPL|^RPS|^MT-", regex=True)
x2 = x2.sort_values('f')

print(
    ggplot(x2)+aes('np.log10(bulk_mat)', 'np.log10(pseudo_mat)')+
        geom_point(aes(color='f'))+
        geom_smooth(method='lm')+
        geom_abline(slope=1, intercept=0)+
        labs(
            x='bulk (RPM)', y='pseudbulkd (RPM)',
            color='RP|MT'
        )
)


# %%
import sparse
import xarray as xa
from combat_atlas.sigs import sigs
from combat_atlas.sigs.fit import fit_gsea

x2 = xa.apply_ufunc(
    lambda x, y: np.corrcoef(x, y)[0,1],
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

#print(
#    ggplot(x2.to_dataframe().reset_index())+aes('R')+geom_freqpoly(bins=30)
#)

x2['t'] = x2.R**2
x2.t.data = sparse.COO(x2.t.data)
x2 = xa.merge([x2, sigs.all1.rename('s')], join='inner')

x5 = fit_gsea(x2.t.expand_dims(g=[1]), x2.s, 1e5)
x5 = x5.squeeze('g').drop('g')
x2.t.data = x2.t.data.todense()

x6 = x5[['ES', 'NES', 'padj']].to_dataframe().reset_index()
x6 = x6.sort_values('padj')

# %%
x6[x6.sig.str.contains('NEUTRO')]
x6.query('NES<0 & padj<0.2')

# %%
x7 = x5.sel(sig='GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_GRAM_NEGATIVE_BACTERIUM')
x7 = x7.sel(gene=x7.leadingEdge.todense()==1)
x7 = x2.sel(gene=x7.gene.data)
x7 = x7.sel(sig=x7.s.sum(dim='gene').todense()>0)

x7[['R', 'symbol']].to_dataframe().sort_values('R')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index()
print(
    ggplot(x2)+aes('clin_Source', 'bulk_mat')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+
        coord_flip()+
        labs(x=f'{x2.feature_id.iloc[0]} bulk (RPM)')+
        theme(legend_position='none')
)

print(
    ggplot(x2)+aes('clin_Source', 'pseudo_mat')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+
        coord_flip()+
        labs(x=f'{x2.feature_id.iloc[0]} psuedobulk (RPM)')+
        theme(legend_position='none')
)

# %%
x3 = self.data1.copy()
x3['bulk_mat'] = 1e6*x3.bulk_mat/x3.bulk_mat.sum(dim='gene_id')
x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index()
print(
    ggplot(x2)+aes('clin_Source', 'pseudo_mat')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+
        facet_grid('cell_type~.', scales='free_y')+
        labs(y=f'{x2.feature_id.iloc[0]} psuedobulk (RPM)')+            
        theme(
            legend_position='none', figure_size=(4, 18*2),
            axis_text_x = element_text(angle=45, hjust=1)
        )
)

# %%
