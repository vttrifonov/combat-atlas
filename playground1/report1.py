# %%
from plotnine import *
import matplotlib.pyplot as plt
import numpy as np
import xarray as xa
from combat_atlas.playground1 import analysis

# %%
self = analysis

# %%
x3 = self.data1.copy()
#x3 = x3.sel(gene_id=(~x3.feature_id.to_series().str.contains("^RPL|^RPS|^MT-", regex=True)).to_xarray())
#x3 = x3.sel(gene_id=~x3.feature_id.isin([
#    'ACTB', 'B2M', 'EEF1A1', 'TMSB4X', 
#    'PTMA','FAU', 'TPT1', 'FTL', 'FTH1','TMSB10', 
#    'HLA-B', 'HLA-A', 'HLA-C'
#]))
#x3['bulk_mat'] = 1e6*x3.bulk_mat/x3.bulk_mat.sum(dim='gene_id')
x3['pseudo_mat'] = x3.pseudo_mat.sum(dim='cell_type')
x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
x3 = x3.drop_dims('cell_type')
x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index() 
print(
    ggplot(x2)+aes('(bulk_mat)', '(pseudo_mat)')+
        geom_point(aes(color='clin_Source'))+
        geom_smooth(method='lm')+
        geom_hline(yintercept=x2.pseudo_mat.mean(), linetype='--')+
        geom_vline(xintercept=x2.bulk_mat.mean(), linetype='--')+
        labs(
            x='bulk (RPM)', y='pseudbulkd (RPM)', 
            title=f'{x2.feature_id.iloc[0]} (R={x2[["bulk_mat", "pseudo_mat"]].corr().to_numpy()[0,1]:.2f})'
        )
)


# %%
x2 = x3.isel(sample_id=1).to_dataframe().reset_index()
x2['f'] = x2.feature_id.str.contains("^RPL|^RPS|^MT-", regex=True)
x2 = x2.sort_values('f')

#x2.query('f==False')[['feature_id', 'bulk_mat', 'pseudo_mat']].sort_values('pseudo_mat').tail(11).feature_id.to_list()

print(
    ggplot(x2)+aes('np.log10(bulk_mat)', 'np.log10(pseudo_mat)')+
        geom_point(aes(color='f'), alpha=0.2)+
        geom_smooth(method='lm')+
        geom_abline(slope=1, intercept=0)+
        labs(
            x='bulk (log10RPM)', y='pseudbulkd (log10RPM)',
            color='RP|MT'
        )
)

# %%
from combat_atlas.sigs import sigs
x5 = self.data3
x6 = x5[['ES', 'NES', 'padj']].to_dataframe().reset_index()
x6 = x6.sort_values('padj')
x6 = x6.query('padj<0.005').reset_index(drop=True)

x6[x6.sig.str.contains('_T_CELL|NEUTRO')].sort_values('ES')

x7 = self.data3.sel(sig='GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION')
x7 = x7.sel(gene=x7.leadingEdge.todense()==1)
#x7 = xa.merge([x7, sigs.all1.rename(sig='sig1')], join='inner')
#x7 = x7.sel(sig1=x7.set.sum(dim='gene').todense()>0)
#x8 = x7.sig1.to_series().reset_index(drop=True)
x7[['symbol', 'R']].to_dataframe().sort_values('R')

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
x2 = x3.sel(gene_id=x3.feature_id=='CD3G').to_dataframe().reset_index()
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
x2 = x3.sel(gene_id=x3.feature_id=='CEACAM8').to_dataframe().reset_index()
print(
    ggplot(x2)+aes('clin_Source', 'np.log10(bulk_mat+1)')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.1)+
        coord_flip()+
        labs(x=f'{x2.feature_id.iloc[0]} bulk (RPM)')+
        theme(legend_position='none')
)

print(
    ggplot(x2)+aes('clin_Source', 'np.log10(pseudo_mat+1)')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+
        coord_flip()+
        labs(x=f'{x2.feature_id.iloc[0]} psuedobulk (RPM)')+
        theme(legend_position='none')
)

# %%
x3 = self.data1.copy()
x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index()
print(
    ggplot(x2)+aes('clin_Source', '(pseudo_mat)')+
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
