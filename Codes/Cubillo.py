import sunpy.map as smap
import glob
import matplotlib.pyplot as plt
import matplotlib

folder = '/home/oan1803/Documents/MOAT/Mancha1/'
files = sorted(glob.glob(folder+'*.fits'))

cubo = smap.Map(files, sequence = True)

cubo.plot()

plt.show()
