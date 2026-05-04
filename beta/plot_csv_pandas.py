import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sys
import os

# ── Config ─────────────────────────────────────────────────────────────────
CSV_FILE = sys.argv[1] if len(sys.argv) > 1 else "data.csv"

# ── Load ────────────────────────────────────────────────────────────────────
if not os.path.exists(CSV_FILE):
    print(f"Error: file '{CSV_FILE}' not found.")
    sys.exit(1)

df = pd.read_csv(CSV_FILE)

# Normalise column names (strip whitespace, lowercase)
df.columns = df.columns.str.strip().str.lower()

if "time" not in df.columns or "width" not in df.columns:
    raise ValueError(f"Expected columns 'time' and 'width', got: {list(df.columns)}")

# Try to parse time as datetime; fall back to numeric
try:
    df["time"] = pd.to_datetime(df["time"])
    time_is_datetime = True
except (ValueError, TypeError):
    df["time"] = pd.to_numeric(df["time"])
    time_is_datetime = False

df = df.sort_values("time")

# ── Plot ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(df["time"], df["width"], linewidth=1.8, color="#2563eb", zorder=3)
ax.fill_between(df["time"], df["width"], alpha=0.12, color="#2563eb")

# Datetime x-axis formatting
if time_is_datetime:
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate()

# Grid & style
ax.grid(True, linestyle="--", alpha=0.5, zorder=0)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)

# Labels
title = os.path.splitext(os.path.basename(CSV_FILE))[0].replace("_", " ").title()
ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
ax.set_xlabel("Time", fontsize=11)
ax.set_ylabel("Value", fontsize=11)

# Stats annotation
stats = (
    f"n={len(df):,}  "
    f"min={df['width'].min():.3g}  "
    f"max={df['width'].max():.3g}  "
    f"mean={df['width'].mean():.3g}"
)
ax.annotate(stats, xy=(0.01, 0.97), xycoords="axes fraction",
            va="top", fontsize=9, color="gray")

plt.tight_layout()

# Save next to the CSV, or in cwd
out_path = os.path.splitext(CSV_FILE)[0] + "_plot.png"
plt.savefig(out_path, dpi=150)
print(f"Saved → {out_path}")
plt.show()
