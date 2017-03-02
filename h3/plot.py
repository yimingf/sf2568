import numpy as np
import matplotlib.pyplot as plt

num_processors = np.array([1., 2., 3., 4., 8., 16.])
elapsed_time = 0.1*np.array([2.850, 1.582, 1.330, 1.108, 1.164, 1.390])

plt.plot(num_processors, elapsed_time)
plt.show()