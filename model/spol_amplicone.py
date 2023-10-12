#%%
#%%
from functools import partial
from random import choices, randint, randrange, random, sample
from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
import pandas as pd
from collections import Counter
from tqdm import tqdm
import time
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from plotly import graph_objects as go
import json
from imp import reload

#%%
spol_list = np.arange(3117003,3127206,1)
weight = [1]*len(spol_list)
spol_data = pd.DataFrame({'spol':spol_list,'weight':weight})
