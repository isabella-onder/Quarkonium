import matplotlib.pyplot as plt

fig = plt.figure()

gs = fig.add_gridspec(
    2, 2,
    width_ratios=[6, 1],
    height_ratios=[1, 1]
)

ax1 = fig.add_subplot(gs[0, 0])   # top-left
ax2 = fig.add_subplot(gs[1, 0])   # bottom-left
ax3 = fig.add_subplot(gs[:, 1])   # right column spanning both rows

plt.subplots_adjust(wspace=0, hspace=0.03)
plt.show()