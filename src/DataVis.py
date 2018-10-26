import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#filename = r"Z:\NE723\FinalProgram\flux_map_A.dat"
#filename = r"Z:\NE723\FinalProgram\flux_map_B.dat"
filename = r"C:\Users\Ryan W\Downloads\SnProgram.dat"


x, y, phi_1 = np.loadtxt(filename).T #Transposed for easier unpacking
nrows, ncols = 10, 10
#nrows, ncols = 151, 151
grid_1 = phi_1.reshape((nrows, ncols))
#grid_2 = phi_2.reshape((nrows, ncols))
#grid_3 = phi_3.reshape((nrows, ncols))
#grid_4 = phi_4.reshape((nrows, ncols))
#grid_5 = phi_5.reshape((nrows, ncols))
#grid_6 = phi_6.reshape((nrows, ncols))
#grid_7 = phi_7.reshape((nrows, ncols))

plt.imshow(grid_1, extent=(x.min(), x.max(), y.max(), y.min()),
           interpolation='nearest', cmap=cm.jet)
plt.colorbar()
plt.savefig('Phi_1.png', bbox_inches='tight')
plt.clf()

# plt.imshow(grid_2, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_2.png', bbox_inches='tight')
# plt.clf()
#
# plt.imshow(grid_3, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_3.png', bbox_inches='tight')
# plt.clf()
#
# plt.imshow(grid_4, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_4.png', bbox_inches='tight')
# plt.clf()
#
# plt.imshow(grid_5, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_5.png', bbox_inches='tight')
# plt.clf()
#
# plt.imshow(grid_6, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_6.png', bbox_inches='tight')
# plt.clf()
#
# plt.imshow(grid_7, extent=(x.min(), x.max(), y.max(), y.min()),
#            interpolation='nearest', cmap=cm.jet)
# plt.colorbar()
# plt.savefig('Phi_7.png', bbox_inches='tight')
# plt.clf()

# f, axarr = plt.subplots(4, 2)
# axarr[0, 0].imshow(grid_1, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[0, 0].set_title('Group 1')
# axarr[0, 1].imshow(grid_2, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[0, 1].set_title('Group 2')
# axarr[1, 0].imshow(grid_3, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[1, 0].set_title('Group 3')
# axarr[1, 1].imshow(grid_4, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[1, 1].set_title('Group 4')
# axarr[2, 0].imshow(grid_5, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[2, 0].set_title('Group 5')
# axarr[2, 1].imshow(grid_6, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[2, 1].set_title('Group 6')
# axarr[3, 0].imshow(grid_7, extent=(x.min(), x.max(), y.max(), y.min()), interpolation='nearest', cmap=cm.jet)
# axarr[3, 0].set_title('Group 7')
# plt.delaxes(axarr[3, 1])
# plt.show()


