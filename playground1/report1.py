# %%
import re
from plotnine import *
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
x2 = x3[['pseudo_mat', 'bulk_mat']].to_dataframe().reset_index()
x2 = x2.groupby(['gene_id', 'feature_id']).apply(
    lambda x: x[['pseudo_mat', 'bulk_mat']].corr().iloc[0,1]
)
x2 = x2.rename('R').reset_index()
print(
    ggplot(x2)+aes('R')+
        geom_freqpoly(30)        
)

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
