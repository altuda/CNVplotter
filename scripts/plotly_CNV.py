import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
import plotly.io as pio
import numpy as np
from scipy import stats



df = pd.read_table('/Users/dasenka/reconCNV/my_samples/patient_07_S152758_R1_merged_val_1_bismark_bt2_pe.deduplicated_sorted.cnr')
chromSizes = pd.read_table('/Users/dasenka/MUW/genome/hg19/hg19.chrom.sizes', header=None, names=['#chrom', 'size']) #chrom sizes

classic_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
chromSizes = chromSizes[chromSizes['#chrom'].isin(classic_chromosomes)]
chromSizes = chromSizes.sort_values('#chrom', key=lambda x: pd.Categorical(x, categories=classic_chromosomes, ordered=True))


print(df.head())
print(chromSizes.head())


# Process chromosome sizes
chromSizes['cumulative_position'] = chromSizes['size'].cumsum() - chromSizes['size']
chromSizes_dict = dict(zip(chromSizes['#chrom'], chromSizes['cumulative_position']))
chromSizes['chrom_center'] = chromSizes['cumulative_position'] + chromSizes['size'] / 2
chrom_centers = chromSizes['chrom_center'].tolist()

# Calculate genomic position
df['genomic_position'] = df['start'] + df['chromosome'].map(chromSizes_dict)

# Filter df to include only classic chromosomes
df = df[df['chromosome'].isin(classic_chromosomes)]

# Calculate chromosome centers
chromSizes['chrom_center'] = chromSizes['cumulative_position'] + chromSizes['size'] / 2
chrom_centers = chromSizes['chrom_center'].tolist()

# Create chromosome boundary lines
chrom_boundaries = chromSizes['cumulative_position'].tolist()[1:] + [chromSizes['size'].sum()]

df['z_log2'] = stats.zscore(df['log2'])

# Set color scale limits based on z-scores
z_min = -3  # Typically, values below -3 are considered outliers
z_max = 3   # Typically, values above 3 are considered outliers

# Create the figure
fig = go.Figure()


# Add scatter points with color gradient
scatter = go.Scattergl(
    x=df['genomic_position'],
    y=df['log2'],  # Keep original log2 values for y-axis
    mode='markers',
    marker=dict(
        color=df['z_log2'],  # Use z-transformed values for color
        colorscale='balance',  # Red-Gray-Blue colorscale
        cmin=z_min,
        cmax=z_max,
        size=3,
        opacity=0.7,
        colorbar=dict(
            title='Z-score of Log2(CNV)',
            thickness=20,
            len=0.5,
            y=0.5,
            yanchor='middle'
        )
    ),
    hovertemplate='<br>'.join([
        'Chromosome: %{customdata[0]}',
        'Start: %{customdata[1]}',
        'Gene: %{customdata[2]}',
        'Depth: %{customdata[3]}',
        'Weight: %{customdata[4]}',
        'Log2: %{y:.3f}',
        'Z-score: %{marker.color:.3f}'
    ]),
    customdata=df[['chromosome', 'start', 'gene', 'depth', 'weight']],
    name='CNV Data'
)
fig.add_trace(scatter)

# Add smoothed line
window_size = 101
df['smooth_log2'] = df['log2'].rolling(window=window_size, center=True).mean()
smoothed_line = go.Scatter(
    x=df['genomic_position'],
    y=df['smooth_log2'],
    mode='lines',
    line=dict(color='purple', width=2),
    name='Smoothed CNV'
)
fig.add_trace(smoothed_line)


y_min = df['log2'].min()
y_max = df['log2'].max()
y_range = [y_min - 0.1 * abs(y_min), y_max + 0.1 * abs(y_max)]  # Add 10% padding
chrom_boundaries = [0] + chromSizes['cumulative_position'].tolist()[1:]

# Customize layout
fig.update_layout(
    title="Genome-Wide Interactive CNV Plot",
    xaxis_title="Chromosome",
    yaxis_title="Log2",
    autosize=True,
    showlegend=False,
    plot_bgcolor='white',
    xaxis=dict(
        tickmode='array',
        tickvals=chrom_centers,  # Use chrom_centers for tick positions
        ticktext=[chrom.replace('chr', '') for chrom in chromSizes['#chrom']],
        tickangle=0,
        tickfont=dict(size=10),
        showgrid=False,  # Remove the grid
        gridcolor='lightgrey',
        gridwidth=1
    ),
    yaxis=dict(
        range=y_range,
        showgrid=True,
        gridcolor='lightgrey',
        gridwidth=0.2
    ),
    margin=dict(l=50, r=50, t=50, b=50),
)

# Add vertical lines to separate chromosomes
for boundary in chrom_boundaries:
    fig.add_shape(
        type="line",
        x0=boundary,
        y0=y_range[0],
        x1=boundary,
        y1=y_range[1],
        line=dict(color="lightgrey", width=1.5, dash='solid'),
    )


# Set the default renderer to 'browser'
pio.renderers.default = "browser"

# Display the plot
fig.show()

print("Plot opened in browser.")

# # Generate the plot content
# plot_div = fig.to_html(full_html=False, include_plotlyjs=False, config={'responsive': True})

# # Create a custom HTML file with proper DOCTYPE and responsive design
# html_content = f"""
# <!DOCTYPE html>
# <html lang="en">
# <head>
#     <meta charset="utf-8">
#     <meta name="viewport" content="width=device-width, initial-scale=1.0">
#     <title>CNV Plot</title>
#     <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
#     <style>
#         body, html {{
#             margin: 0;
#             padding: 0;
#             height: 100%;
#             width: 100%;
#         }}
#         #plotly-div {{
#             width: 100%;
#             height: 100vh;
#         }}
#     </style>
# </head>
# <body>
#     <div id="plotly-div">{plot_div}</div>
#     <script>
#         window.onload = function() {{
#             var gd = document.getElementById('plotly-div');
#             Plotly.relayout(gd, {{
#                 'xaxis.autorange': true,
#                 'yaxis.autorange': true
#             }});
#         }};
#         window.addEventListener('resize', function() {{
#             Plotly.Plots.resize(document.getElementById('plotly-div'));
#         }});
#     </script>
# </body>
# </html>
# """

# # Write the HTML content to a file
# with open("full_screen_plot.html", "w") as f:
#     f.write(html_content)

# print("HTML file generated successfully.")