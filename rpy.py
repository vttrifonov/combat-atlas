# %%
import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects import r as R
import xarray as xa
import numpy as np
import pandas as pd

dict2ri = ro.conversion.Converter('dict converter')
@dict2ri.py2rpy.register(dict)    
def _dict_py2rpy(x):
    return ro.ListVector(x)

@dict2ri.rpy2py.register(ri.ListSexpVector)
def _dict_rpy2py(x):
    co = ro.conversion.get_conversion()
    with localconverter(ro.default_converter) as co1:
        x = co1.rpy2py(x)
        return {
            k: co.rpy2py(co1.rpy2py(v)) 
            for k, v in x.items()
        }
    
xa2ri = ro.conversion.Converter('xarray converter')

@xa2ri.py2rpy.register(xa.DataArray)   
def _xa_py2rpy(d):
    with localconverter(ro.default_converter+numpy2ri.converter+dict2ri) as co:
        array = co.py2rpy(np.asarray(d.data, order='C'))
        dimnames = co.py2rpy({k: d[k].data for k in d.dims})
        dims = co.py2rpy(np.asarray(d.data.shape))
    with localconverter(ro.default_converter):
        array = R.array(array, dim=dims)
        array.dimnames = dimnames
    return array

@xa2ri.rpy2py.register(ro.Array)
def _xa_rpy2py(x):
    with localconverter(ro.default_converter+numpy2ri.converter) as co:
        dimnames = co.rpy2py(x.dimnames)
        matrix = co.rpy2py(x)
    return xa.DataArray(matrix, coords=dimnames)

num2ri = ro.conversion.Converter('num converter')

@num2ri.py2rpy.register(object)
def _num_py2rpy(x):
    conv = ro.default_converter
    if isinstance(x, pd.DataFrame) or isinstance(x, pd.Series):
        conv = conv+pandas2ri.converter
    elif isinstance(x, xa.DataArray):
        conv = conv+xa2ri
    elif isinstance(x, np.ndarray):
        conv = conv+numpy2ri.converter
    with localconverter(conv) as co:
        return co.py2rpy(x)

@num2ri.rpy2py.register(object)
def _num_rpy2py(x):
    with localconverter(xa2ri+dict2ri) as co:        
        return co.rpy2py(x)
    
# %%
def r_func(f):
    from rpy2.robjects.conversion import localconverter
    def wrap(*args):
        with localconverter(num2ri+dict2ri):
            return f(*args)
    return wrap