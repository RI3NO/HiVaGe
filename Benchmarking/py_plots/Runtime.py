import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from CSV files into DataFrames
Mair_df = pd.read_csv('Data/Mair_v3_13082023_metrics.csv')
Richard_df = pd.read_csv('Data/Richard_v2_13082023_metrics (2).csv')
Campbell_df = pd.read_csv('Data/Campbell_v1_13082023_metrics.csv')
Buettner_df = pd.read_csv('Data/Buettner_v3_14082023_metrics.csv')
HIV_df = pd.read_csv('Data/HIV_v1_13082023_metrics.csv')

# Select specific columns of interest
col_of_interest = ['Unnamed: 0', 'Runtime']
Mair_df = Mair_df[col_of_interest]
Richard_df = Richard_df[col_of_interest]
Campbell_df = Campbell_df[col_of_interest]
Buettner_df = Buettner_df[col_of_interest]
HIV_df = HIV_df[col_of_interest]

# Merge DataFrames
df = Mair_df.merge(Richard_df, on='Unnamed: 0', how='outer')
df = df.merge(Campbell_df, on='Unnamed: 0', how='outer')
df = df.merge(Buettner_df, on='Unnamed: 0', how='outer')
df = df.merge(HIV_df, on='Unnamed: 0', how='outer')

# Sort columns by 'Method'
df.columns = ['Method', 'Mair', 'Richard', 'Campbell', 'Buettner', 'HIV']
df = df.sort_values(by='Method')
print(df)


# HEATMAP1
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator

plt.figure(figsize=(6.2, 7))
heatmap = sns.heatmap(df.set_index('Method'), cmap="OrRd", annot=False,
                      fmt=".2f", norm=LogNorm(), linewidths=0.5)

cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(which='minor', length=0)


tick_labels = ["NULL", "NULL", '1', '10', r'$10^2$', r'$10^3$']
cbar.ax.set_yticklabels(tick_labels)
heatmap.set_ylabel("")

filename_svg = "Runtime/Runtime_new_heatmap.svg"
plt.savefig(filename_svg, format='svg', bbox_inches='tight')
plt.show()


# LINEPLOT
df_long = pd.melt(df, id_vars=['Method'])
value = np.log(df_long['value'])


ax1 = sns.lineplot(data=df_long, x='Method', y=value, hue='variable')
ax1.set_xticklabels([])
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.legend(title="Dataset")

# Split and save values into an Excel file
list_values = value.to_list()
split_values = np.array_split(list_values, 5)

df_split = pd.DataFrame(split_values, index=[f"Part {i+1}" for i in range(5)])
excel_filename = "split_values.xlsx"
df_split.to_excel(excel_filename, index=False)

print(df['Method'].to_list())

filename_svg = "Runtime/Runtime_lineplot.svg"
plt.savefig(filename_svg, format='svg', bbox_inches='tight')
plt.show()


# HEATMAP2
# Load data from an Excel file for the second heatmap
excel_file_path = "heatmap2.xlsx"
df = pd.read_excel(excel_file_path, index_col=0)
color = 'PuBu'

plt.figure(figsize=(5, 7))
heatmap = sns.heatmap(df, cmap=color, annot=False, fmt=".2f", linewidths=0.5, norm=LogNorm())
tick_labels = ["NULL", "NULL", '1', '10', r'$10^2$', r'$10^3$']
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(which='minor', length=0)
cbar.ax.set_yticklabels(tick_labels)

filename_svg = "Runtime/Runtime_to_Mair.svg"
plt.savefig(filename_svg, format='svg', bbox_inches='tight')
plt.show()