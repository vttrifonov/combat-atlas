#%%
if __name__ == '__main__':
    __package__ = 'ards_geo'

from pathlib import Path
import xarray as xa
import numpy as np
import pandas as pd
import itertools as it

#%%
class _config:
    project='combat-atlas'
    cache = Path.home()/'.cache'/project
    root = Path(__file__).parent

config = _config()

#%%