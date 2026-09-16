#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8.5), dpi=300)

# --- Mode 1: Genomic Position Mode (-w 10000 -l 2000) ---
ax1.set_title("1. Genomic Position Mode (-w 10000 -l 2000)", fontsize=13, fontweight='bold', pad=12)
ax1.set_xlim(0, 16000)
ax1.set_ylim(0, 4.2)
ax1.axis('off')

# Chromosome axis
ax1.annotate('', xy=(15000, 3.2), xytext=(0, 3.2),
             arrowprops=dict(arrowstyle="->", lw=2, color='black'))
ax1.text(7500, 3.75, "Chromosome Base Pairs (bp)", ha='center', fontsize=10, fontweight='bold')

# Ticks
for pos in [1, 2001, 4001, 10000, 12000, 14000]:
    ax1.plot([pos, pos], [3.1, 3.3], color='black', lw=1.5)
    ax1.text(pos, 3.45, f"{pos}", ha='center', fontsize=8)

# Variant positions on axis
vars_dense1 = [200, 400, 600, 900, 1100, 1300, 1600, 1800, 2100, 2400, 2700, 3100, 3400, 3700]
vars_sparse = [4500, 6200, 8100, 9500]
vars_dense2 = [10200, 10600, 11100, 11500, 12100, 12500, 13000, 13400, 13900, 14300, 14700]
all_vars = vars_dense1 + vars_sparse + vars_dense2

for vpos in all_vars:
    ax1.plot([vpos, vpos], [3.15, 3.25], color='red', lw=1)

# Windows in Mode 1 (Fixed 10,000 bp width, variable variant counts)
# Window 1
rect1 = patches.Rectangle((1, 2.1), 9999, 0.6, linewidth=1.5, edgecolor='#2b5c8f', facecolor='#6baed6', alpha=0.7)
ax1.add_patch(rect1)
ax1.text(5000, 2.4, "Window 1: bp 1 - 10000  [Fixed 10,000 bp width  |  18 variants]", ha='center', va='center', fontsize=9, color='white', fontweight='bold')

# Window 2
rect2 = patches.Rectangle((2001, 1.2), 9999, 0.6, linewidth=1.5, edgecolor='#d94701', facecolor='#fdae6b', alpha=0.7)
ax1.add_patch(rect2)
ax1.text(7000, 1.5, "Window 2: bp 2001 - 12000  [Fixed 10,000 bp width  |  12 variants]", ha='center', va='center', fontsize=9, color='black', fontweight='bold')

# Window 3
rect3 = patches.Rectangle((4001, 0.3), 9999, 0.6, linewidth=1.5, edgecolor='#238b45', facecolor='#74c476', alpha=0.7)
ax1.add_patch(rect3)
ax1.text(9000, 0.6, "Window 3: bp 4001 - 14000  [Fixed 10,000 bp width  |  13 variants]", ha='center', va='center', fontsize=9, color='white', fontweight='bold')


# --- Mode 2: Fixed Variant Width Mode (-w 10000 -l 2000 -f) ---
ax2.set_title("2. Fixed Variant Width Mode (-w 10000 -l 2000 -f)", fontsize=13, fontweight='bold', pad=12)
ax2.set_xlim(0, 16000)
ax2.set_ylim(0, 4.2)
ax2.axis('off')

# Chromosome axis
ax2.annotate('', xy=(15000, 3.2), xytext=(0, 3.2),
             arrowprops=dict(arrowstyle="->", lw=2, color='black'))
ax2.text(7500, 3.75, "Chromosome Variants", ha='center', fontsize=10, fontweight='bold')

# Variant positions on axis (well spaced so text doesn't overlap)
ax2.text(100, 3.42, "v_1\n(pos 100)", ha='center', fontsize=8)
ax2.text(1200, 3.42, "v_2001\n(pos 1200)", ha='center', fontsize=8)
ax2.text(2500, 3.42, "v_4001\n(pos 2500)", ha='center', fontsize=8)
ax2.text(4500, 3.42, "v_10000\n(pos 4500)", ha='center', fontsize=8)
ax2.text(12500, 3.42, "v_12000\n(pos 12500)", ha='center', fontsize=8)
ax2.text(15500, 3.42, "v_14000\n(pos 15500)", ha='center', fontsize=8)

for pos in [100, 1200, 2500, 4500, 12500, 15500]:
    ax2.plot([pos, pos], [3.1, 3.3], color='red', lw=1.5)

# Windows in Mode 2 (Fixed 10,000 variants, variable physical bp width)
# Window 1: var 1..10000 (dense region: pos 100 - 4500 -> 4,400 bp wide)
rect1_f = patches.Rectangle((100, 2.1), 4400, 0.6, linewidth=1.5, edgecolor='#2b5c8f', facecolor='#6baed6', alpha=0.7)
ax2.add_patch(rect1_f)
ax2.text(2300, 2.4, "Window 1: 10,000 variants (4,400 bp wide: pos 100 - 4500)", ha='center', va='center', fontsize=8, color='white', fontweight='bold')

# Window 2: var 2001..12000 (sparse region: pos 1200 - 12500 -> 11,300 bp wide)
rect2_f = patches.Rectangle((1200, 1.2), 11300, 0.6, linewidth=1.5, edgecolor='#d94701', facecolor='#fdae6b', alpha=0.7)
ax2.add_patch(rect2_f)
ax2.text(6850, 1.5, "Window 2: 10,000 variants (11,300 bp wide: pos 1200 - 12500)", ha='center', va='center', fontsize=9, color='black', fontweight='bold')

# Window 3: var 4001..14000 (very sparse region: pos 2500 - 15500 -> 13,000 bp wide)
rect3_f = patches.Rectangle((2500, 0.3), 13000, 0.6, linewidth=1.5, edgecolor='#238b45', facecolor='#74c476', alpha=0.7)
ax2.add_patch(rect3_f)
ax2.text(9000, 0.6, "Window 3: 10,000 variants (13,000 bp wide: pos 2500 - 15500)", ha='center', va='center', fontsize=9, color='white', fontweight='bold')

plt.tight_layout()
plt.savefig("window_schematic.png", dpi=300)
print("Updated window_schematic.png successfully")
