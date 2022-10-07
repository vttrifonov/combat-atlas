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
            title=f'{x2.feature_id.iloc[0]} (R={x2[["bulk_mat", "pseudo_mat"]].corr().to_numpy()[0,1]:.2f})',
            color='Source'
        )
)


# %%
x2 = x3.isel(sample_id=1).to_dataframe().reset_index()
x2['f'] = np.where(
    x2.feature_id.str.contains("^RPL|^RPS", regex=True), 'RP',
    np.where(
        x2.feature_id.str.contains("^MT-", regex=True), 'MT',
        np.where(
            x2.feature_id.isin(['ACTB', 'B2M']), x2.feature_id,
            np.nan    
        )
    )
)

x2 = x2.sort_values('f', na_position='first')

#x2.query('f==False')[['feature_id', 'bulk_mat', 'pseudo_mat']].sort_values('pseudo_mat').tail(11).feature_id.to_list()

print(
    ggplot(x2)+aes('bulk_mat', 'pseudo_mat')+
        geom_point(aes(color='f'), alpha=0.8)+
        geom_smooth(method='lm')+
        geom_abline(slope=1, intercept=0)+
        labs(
            x='bulk (RPM)', y='pseudbulkd (RPM)',
            color='mRNA'
        )+
        scale_x_log10()+scale_y_log10()
)

# %%
def _():
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
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1')
x2 = x2.to_dataframe().reset_index()
x2 = x2[['sample_id', 'feature_id', 'pseudo_mat', 'bulk_mat', 'clin_Source']]
x2 = x2.melt(id_vars=['sample_id', 'feature_id', 'clin_Source'])
x2['variable'] = np.where(x2.variable=='pseudo_mat', 'pseudobulk', 'bulk')

print(
    ggplot(x2)+aes('clin_Source', 'value')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+        
        facet_grid('~variable', scales='free_x')+
        labs(y=f'{x2.feature_id.iloc[0]} (RPM)', x='Source')+
        coord_flip()+
        theme(legend_position='none')
)

# %%
x2 = x3.sel(gene_id=x3.feature_id=='CD3G')
x2 = x2.to_dataframe().reset_index()
x2 = x2[['sample_id', 'feature_id', 'pseudo_mat', 'bulk_mat', 'clin_Source']]
x2 = x2.melt(id_vars=['sample_id', 'feature_id', 'clin_Source'])
x2['variable'] = np.where(x2.variable=='pseudo_mat', 'pseudobulk', 'bulk')

print(
    ggplot(x2)+aes('clin_Source', 'value')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+        
        facet_grid('~variable', scales='free_x')+
        labs(y=f'{x2.feature_id.iloc[0]} (RPM)', x='Source')+
        coord_flip()+
        theme(legend_position='none')+
        scale_y_log10()
)

# %%
x2 = x3.sel(gene_id=x3.feature_id=='CEACAM8')
x2 = x2.to_dataframe().reset_index()
x2 = x2[['sample_id', 'feature_id', 'pseudo_mat', 'bulk_mat', 'clin_Source']]
x2 = x2.melt(id_vars=['sample_id', 'feature_id', 'clin_Source'])
x2['variable'] = np.where(x2.variable=='pseudo_mat', 'pseudobulk', 'bulk')

print(
    ggplot(x2)+aes('clin_Source', 'value')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+        
        facet_grid('~variable', scales='free_x')+
        labs(y=f'{x2.feature_id.iloc[0]} (log10RPM)', x='Source')+
        coord_flip()+
        theme(legend_position='none')+
        scale_y_log10()
)


# %%
x3 = self.data1.copy()
x3['pseudo_mat'] = 1e6*x3.pseudo_mat/x3.pseudo_mat.sum(dim='gene_id')
x3 = x3.sel(sample_id=x3.clin_Source!='COVID_HCW_MILD')

# %%
x2 = x3.sel(gene_id=x3.feature_id=='C5AR1').to_dataframe().reset_index()
x2['cell_type'] = np.where(x2.cell_type.str.contains('Mono'), x2.cell_type, 'Other')
print(
    ggplot(x2)+aes('clin_Source', '(pseudo_mat)')+
        geom_violin(aes(fill='clin_Source'))+
        geom_jitter(width=0.2)+
        facet_grid('cell_type~.', scales='free_y')+
        labs(y=f'{x2.feature_id.iloc[0]} psuedobulk (RPM)')+            
        theme(
            legend_position='none', figure_size=(10, 5),
            axis_text_x = element_text(angle=45, hjust=1)
        )+
        scale_y_log10()        
)

# %%
