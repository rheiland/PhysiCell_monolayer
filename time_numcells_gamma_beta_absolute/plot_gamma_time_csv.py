import csv
import matplotlib.pyplot as plt

csv_file = "gamma_time_10K.csv"

bvals = []
tvals = []
critical_times = [15, 27, 68, 543]
print(f"  critical_times = {critical_times}")
with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        bvals.append(float(row['gamma']) * 100)
        tvals.append(float(row['time']))
        num = len(tvals)
        if num > 1:
            for ct in critical_times:
                t_prev = tvals[num-2]
                t = tvals[num-1]
                if ct > t_prev and ct < t:
                    # print(f"len(tvals)= {num}, last t={tvals[-1]}")
                    # print(f"t_pre={tvals[num-2]}, t={tvals[num-1]}")
                    print(f"{t_prev} < {ct} < {t}")

        # for time in critical_times:
            # if time < tvals[-1] 

# tvals *= 100
fig, ax0 = plt.subplots()
# ax0.plot(bvals, tvals, marker='.', markersize=3,color='k')
ax0.plot(bvals, tvals, color='k')

ax0.set_xlim(left=0, right=100)
ax0.set_ylim(bottom=0)
ax0.set_title(r"$\beta=0$", fontsize=12)
# ax0.set_xlabel(r'$\gamma$ (% of $\mathit{f}$-distribution at $10^3$ cells)')
ax0.set_xlabel(r'$\gamma$')
ax0.set_ylabel('Time (calibrated for 5T)')
plt.tight_layout()
plt.savefig('physicell_gamma_time_10K.png')
plt.show()
