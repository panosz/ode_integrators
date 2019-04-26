import pandas as pd
import numpy as np


n_extra_columns = 8
foo_table = pd.read_csv('my_foo.txt',sep=' ')

data = foo_table.values

extra_columns_shape = (data.shape[0],n_extra_columns)

extra_columns = np.zeros(extra_columns_shape)

new_data = np.append(data,extra_columns,axis=1)

np.savetxt('my_foo_extended.txt',new_data,delimiter=' ')
