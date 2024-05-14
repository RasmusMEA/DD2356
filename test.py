import struct
import numpy as np
import matplotlib.pyplot as plt

def test():
	fig = plt.figure(figsize=(4,4), dpi=80)
	f = open("output.bin", "rb")
	# num_doubles = 5 * 5
	dim = 16
	for _ in range(100):
		mtx_list = []
		for _ in range(dim):
			row = f.read(dim * 8)
			format_string = f"{dim}d"
			doubles_row= list(struct.unpack(format_string, row))
			mtx_list.append(doubles_row)

		mtx = np.stack(mtx_list, axis=0)

		plt.cla()
		plt.imshow(mtx.T)
		plt.clim(0.8, 2.2)
		ax = plt.gca()
		ax.invert_yaxis()
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)	
		ax.set_aspect('equal')	
		plt.pause(0.0001)


		


	


test()
