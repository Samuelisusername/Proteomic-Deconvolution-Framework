import matplotlib.pyplot as plt
import numpy as np

# Data
labels = ["1,2 vs 3,4", "3,4 vs 2,1", "2,3 vs 1,4", "1,3 vs 2,4", "4,2 vs 1,3", "1,4 vs 2,3"]
errors = np.array([
    [1.6, 3.3, 3.31],
    [5.5, 6.5, 7.44],
    [4.3, 7.44, 3.1],
    [2.6, 3.09, 7.39],
    [4.34, 3.3, 6.6],
    [2.63, 3.15, 6.66]
])

# Per-bar custom labels
bar_labels = [
    ["1,2→sig\n3,4→bulk", "1→sig\n3,4→bulk", "2→sig\n3,4→bulk"],
    ["3,4→sig\n2,1→bulk", "4→sig\n2,1→bulk", "3→sig\n2,1→bulk"],
    ["2,3→sig\n1,4→bulk", "2→sig\n1,4→bulk", "3→sig\n1,4→bulk"],
    ["1,3→sig\n2,4→bulk", "1→sig\n2,4→bulk", "3→sig\n2,4→bulk"],
    ["4,2→sig\n1,3→bulk", "2→sig\n1,3→bulk", "4→sig\n1,3→bulk"],
    ["1,4→sig\n2,3→bulk", "1→sig\n2,3→bulk", "4→sig\n2,3→bulk"]
]

# Plot
x = np.arange(len(labels))
width = 0.25
fig, ax = plt.subplots(figsize=(14, 8))

# Bar containers
bars1 = ax.bar(x - width, errors[:, 0], width, color='#1f77b4', label='Two Reference \nPatients  for \nSignature Matrix')
bars2 = ax.bar(x,         errors[:, 1], width, color='gray', label='One Reference \nPatient for \nSignature Matrix ')
bars3 = ax.bar(x + width, errors[:, 2], width, color='gray', label='Other Reference\n Patient for \nSignature Matrix')

# Annotate each bar with height and label underneath
for i in range(len(x)):
    for j, bars in enumerate([bars1, bars2, bars3]):
        bar = bars[i]
        height = bar.get_height()

        # Add value label on top
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=12)

        # Add per-bar label under bar
        ax.annotate(bar_labels[i][j],
                    xy=(bar.get_x() + bar.get_width()/2, 0),
                    xytext=(0, -25), textcoords="offset points",
                    ha='center', va='top', fontsize=10,
                    rotation=45)

# Style
ax.set_ylabel('Mean Absolute Error [%]', fontsize=20)
ax.set_title('Two vs. One Signature Matrix Performance (MAE)', fontsize = 18)
ax.set_xticks(x)
#ax.set_xticklabels(labels, rotation=45)
ax.set_ylim(0, max(errors.flatten()) + 1)
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=15)

# Add margin to bottom to fit rotated labels
plt.subplots_adjust(bottom=0.3)
plt.tight_layout()
plt.show()