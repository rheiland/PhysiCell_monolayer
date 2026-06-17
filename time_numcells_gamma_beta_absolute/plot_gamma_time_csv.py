import csv
import matplotlib.pyplot as plt

csv_file = "gamma_time_10K.csv"

bvals = []
tvals = []
with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        bvals.append(float(row['gamma']) * 100)
        tvals.append(float(row['time']))

# tvals *= 100
fig, ax0 = plt.subplots()
# ax0.plot(bvals, tvals, marker='.', markersize=3,color='k')
ax0.plot(bvals, tvals, color='k')

ax0.set_xlim(left=0, right=100)
ax0.set_ylim(bottom=0)
ax0.set_title(r"$\beta=0$", fontsize=12)
ax0.set_xlabel(r'$\gamma$ (% of $\mathit{f}$-distribution at $10^3$ cells)')
ax0.set_ylabel('Time (calibrated for 5T)')
plt.tight_layout()
plt.savefig('physicell_gamma_time_10K.png')
plt.show()
