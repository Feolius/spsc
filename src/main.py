import numpy as np

a = np.array([1, 2, 3, 4, 5, 6])

b = a.reshape((2,3))
b = b.reshape(a.size)
c = a.shape
print b[2]
