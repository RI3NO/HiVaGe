import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from CSV files
Mair_df = pd.read_csv('Data/Mair_v3_13082023_metrics.csv')
Richard_df = pd.read_csv('Data/Richard_v2_13082023_metrics (2).csv')
Campbell_df = pd.read_csv('Data/Campbell_v1_13082023_metrics.csv')
Buettner_df = pd.read_csv('Data/Buettner_v3_14082023_metrics.csv')
HIV_df = pd.read_csv('Data/HIV_v1_13082023_metrics.csv')

# Store the DataFrames in a list
datasets = [Mair_df, Richard_df, Campbell_df, Buettner_df, HIV_df]

fig, axes = plt.subplots(nrows=1, ncols=len(datasets), figsize=(20, 3))

# Iterate through each subplot, dataset, and its name
for ax, dataset_name, dataset in zip(axes, ['Mair', 'Richard', 'Campbell', 'Buettner', 'HIV'], datasets):

    # Extract specific columns from the dataset
    calinski = dataset['Calinski-Harabasz index']
    over_h = dataset['Overlap with highly expressed genes'].to_list()
    over_l = dataset['Overlap with lowly expressed genes'].to_list()

    # Define colors
    back_color = "#C3CFF8"
    line_color = "#1C008D"

    # Create a bar plot and line plots
    sns.barplot(y=over_h, x=np.arange(0, len(calinski), step=1), color=back_color, ax=ax)
    sns.lineplot(y=over_l, x=np.arange(0, len(calinski), step=1), color=line_color, ax=ax)
    ax2 = ax.twinx()
    sns.lineplot(calinski, color="#1082D4", ax=ax2)

    # Customize subplot appearance
    ax.set_xticklabels([])
    ax.set_ylabel("")
    ax2.set_ylabel("")
    ax.set_title(dataset_name)

# Save the plots
plt.tight_layout()
filename_svg = "cal_over/all_datasets_Calinski_OverH.svg"
filename_jpg = "cal_over/all_datasets_Calinski_OverH.jpg"
plt.savefig(filename_svg, format='svg')
plt.savefig(filename_jpg, format='jpg')
plt.show()
