import csv
import matplotlib.pyplot as plt

csv_file = "beta_time_10K.csv"

bvals = []
tvals = []
with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # bvals.append(float(row['beta']) * 100)
        bvals.append(float(row['beta']))
        tvals.append(float(row['time']))

fig, ax0 = plt.subplots()
# ax0.plot(bvals, tvals, marker='.', markersize=3,color='k')
# ax0.plot(bvals, tvals, color='k')
ax0.plot(bvals, tvals, '.-', color='lightskyblue', marker='o', markersize=4)

# ax0.set_xlim(left=0, right=100)
ax0.set_xlim(left=-0.05, right=1.05)
ax0.set_title(r"$\beta$ sweep ($\gamma=0$)", fontsize=12)
# ax0.set_xlabel(r'$\beta$ (% of $a$-distribution at $10^3$ cells)')
#ax0.set_xlabel(r'$\beta$ (% of $\mathit{a}$-distribution at $10^3$ cells)')
ax0.set_xlabel(r'$\beta$')
ax0.set_yscale('log')
ax0.set_ylim(0, 10**3) 

ax0.set_ylabel(r"Time to reach $10^4$ cells - (Cell Cycle Time units)")
plt.tight_layout()
plt.savefig('physicell_beta_time_10K.png')
plt.show()
