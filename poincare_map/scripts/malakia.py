import numpy as np
from functools import partial


#  @partial(np.vectorize, excluded=['b'])
@np.vectorize
def myfunc(a, b):
    if a > b:
        return a-b
    else:
        return a+b


print(myfunc([[2, 3], [4, 5]], [1, 3]))
