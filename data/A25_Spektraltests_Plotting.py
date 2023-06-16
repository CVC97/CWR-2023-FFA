import matplotlib.pyplot as plt
import numpy as np


# Aufgabe 25 - Spektraltests Zufallsgeneratoren
A25_Randu_data = np.loadtxt("./A25_Randu.csv", delimiter = ",", skiprows = 0)
A25_Mersenne_Twister_data = np.loadtxt("./A25_Mersenne_Twister.csv", delimiter = ",", skiprows = 0)

randu_x = A25_Randu_data[:,0]
randu_y = A25_Randu_data[:,1]
randu_z = A25_Randu_data[:,2]

mt_x = A25_Mersenne_Twister_data[:,0]
mt_y = A25_Mersenne_Twister_data[:,1]
mt_z = A25_Mersenne_Twister_data[:,2]


fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection = '3d')
ax.plot(randu_x, randu_y, randu_z, '.', color = 'xkcd:red pink', alpha = 0.25, label = 'RANDU')
ax.plot(mt_x, mt_y, mt_z, '.', color = 'blue', alpha = 0.25, label = 'Mersenne Twister')

ax.legend()
ax.legend(loc="upper right")

plt.show()