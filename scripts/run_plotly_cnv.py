import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import plotly.express as px
from plotly.subplots import make_subplots
from scipy import stats
from typing import List, Dict, Tuple
import argparse
import os
import logging

# Constants
CLASSIC_CHROMOSOMES = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
DEFAULT_GENES = [
    "MDM4", "MYCN", "IDH1", "GLI2", "PIK3CA", "CTNNB1", "PDGFRA", "FGFR3", "TERT", "APC",
    "MYB", "BRAF", "EGFR", "CDK6", "MET", "MYC", "MYBL1", "CDKN2A", "CDKN2B",
    "PTCH1", "PTEN", "SUFU", "MGMT", "RELA", "CCND1", "CCND2", "CDK4", "MDM2",
    "RB1", "BRCA2", "IDH2", "CREBBP", "TSC2", "TP53", "NF1", "PPM1D",
    "SMARCA4", "GNAS", "NF2", "SMARCB1"
]
SIZE_MIN = 2
SIZE_MAX = 10
SMOOTHING_WINDOW_SIZE = 101  # Window size for rolling mean

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate CNV plots for multiple samples.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing CNV files')
    parser.add_argument('-c', '--chrom_sizes', required=True, help='Path to chromosome sizes file')
    parser.add_argument('-o', '--output_dir', default='output', help='Directory for output files')
    parser.add_argument('-p', '--prefix', default='CNV_plot', help='Prefix for output file names')
    parser.add_argument('-g', '--genes', help='Comma-separated list of additional genes to highlight')
    return parser.parse_args()

def read_data(file_path: str, is_chrom_sizes: bool = False, raise_errors: bool = True) -> pd.DataFrame:
    """
    Read data from a file, optionally processing chromosome sizes, with configurable error handling.

    Args:
        file_path (str): Path to the data file.
        is_chrom_sizes (bool): Flag to indicate if the file contains chromosome sizes.
        raise_errors (bool): Flag to indicate if errors should be raised or logged.

    Returns:
        pd.DataFrame: The loaded data.
    """
    try:
        # Read the data file
        df = pd.read_table(
            file_path,
            names=['#chrom', 'size'] if is_chrom_sizes else None,
            sep='\s+' if is_chrom_sizes else None,
            engine='python'
        )
        if df.empty:
            error_msg = f"No data found in {file_path}"
            if raise_errors:
                raise ValueError(error_msg)
            else:
                logging.error(error_msg)
                return pd.DataFrame()  # Return empty DataFrame

        if is_chrom_sizes:
            logging.info("Processing chromosome sizes in read_data")
            df = process_chrom_sizes(df)

        return df

    except Exception as e:
        error_msg = f"Error reading {file_path}: {e}"
        if raise_errors:
            logging.error(error_msg, exc_info=True)
            raise RuntimeError(f"Failed to read {file_path}") from e
        else:
            logging.error(error_msg)
            return pd.DataFrame()  # Return empty DataFrame

def process_chrom_sizes(df: pd.DataFrame) -> pd.DataFrame:
    """
    'cumulative_position' - calculate cumulative sum, shift(1) then shifts the values down by 1
    (since the cumulative position of a chromosome should start at the sum of the sizes of all previous chromosomes).
    The fillna(0) replaces the NaN in the first row (resulting from the shift) with 0.
    """
    df = df[df['#chrom'].isin(CLASSIC_CHROMOSOMES)]
    df = df.sort_values('#chrom', key=lambda x: pd.Categorical(x, categories=CLASSIC_CHROMOSOMES, ordered=True))
    df = df.reset_index(drop=True)
    df['cumulative_position'] = df['size'].cumsum().shift(1).fillna(0)
    df['chrom_center'] = df['cumulative_position'] + df['size'] / 2
    return df

def get_highlight_genes(args_genes: str) -> Tuple[List[str], List[str]]:
    default_genes = DEFAULT_GENES.copy()
    custom_genes = [gene.strip() for gene in args_genes.split(',')] if args_genes else []
    return default_genes, custom_genes

def calculate_genomic_position(df: pd.DataFrame, chrom_sizes: pd.DataFrame) -> pd.DataFrame:
    """
    1. zipping two columns from chrom_sizes and create a dictionary
    2. create genomic position column in df by adding the start position of the data point to the cumulative position of the chromosome.
    The cumulative position is retrieved from the chromSizes_dict using the chromosome name as the key (get()).
    The axis=1 parameter specifies that the function should be applied to rows, not columns.
    """
    chromSizes_dict = dict(zip(chrom_sizes['#chrom'], chrom_sizes['cumulative_position']))
    df['genomic_position'] = df.apply(lambda row: row['start'] + chromSizes_dict.get(row['chromosome'], 0), axis=1)
    return df

def create_scatter_plot(
    df: pd.DataFrame,
    default_genes: List[str],
    custom_genes: List[str]
) -> Tuple[List[go.Scattergl], Dict[str, List[Dict]]]:
    """
    Creates scatter plot traces and prepares annotations for default and custom genes.

    Returns:
        traces: List of Plotly Scattergl traces.
        annotations: Dictionary containing separate annotations for default and custom genes.
    """
    # Calculate z-scores and point sizes
    df['z_log2'] = stats.zscore(df['log2'])
    df['point_size'] = SIZE_MIN + df['weight'] * (SIZE_MAX - SIZE_MIN)

    # Main scatter trace
    main_scatter = go.Scattergl(
        x=df['genomic_position'],
        y=df['log2'],
        mode='markers',
        marker=dict(
            color=df['z_log2'],
            colorscale='balance',
            cmin=-3, #determined by z-score meaning, 3 is really an outlier
            cmax=3,
            size=df['point_size'],
            sizemode='diameter',
            opacity=0.7,
            line=dict(width=0),
            colorbar=dict(
                title='Z-score of Log2 ratio',
                thickness=20,
                lenmode='fraction',
                len=0.5,
                y=0.85,
                yanchor='top'
            ),
        ),
        hovertemplate='<br>'.join([
            'Chromosome: %{customdata[0]}',
            'Start: %{customdata[1]}',
            'Gene: %{customdata[2]}',
            'Depth: %{customdata[3]}',
            'Weight: %{customdata[4]:.3f}',
            'Log2: %{y:.3f}',
            'Z-score: %{marker.color:.3f}'
        ]),
        customdata=df[['chromosome', 'start', 'gene', 'depth', 'weight']],
        name='CNV Data'
    )

    # Smoothed line trace
    df['smooth_log2'] = df['log2'].rolling(window=SMOOTHING_WINDOW_SIZE, center=True).mean()
    smoothed_line = go.Scattergl(
        x=df['genomic_position'],
        y=df['smooth_log2'],
        mode='lines',
        line=dict(color='purple', width= 2),
        name='Smoothed',
        hoverinfo='skip'
    )

    traces = [main_scatter, smoothed_line]
    annotations = {'default': [], 'custom': []}

    for gene_type, gene_list in [('default', default_genes), ('custom', custom_genes)]:
        highlighted_df = df[df['gene'].isin(gene_list)].copy()
        if not highlighted_df.empty:
            color = 'red' if gene_type == 'default' else '#246A73'
            traces.append(create_highlighted_scatter(highlighted_df, color, f"Highlighted {gene_type.capitalize()} Genes"))
            annotations[gene_type] = create_gene_annotations(highlighted_df, color)

    return traces, annotations

def create_highlighted_scatter(df: pd.DataFrame, color: str, name: str) -> go.Scattergl:
    return go.Scattergl(
        x=df['genomic_position'],
        y=df['log2'],
        mode='markers',
        marker=dict(
            color=color,
            size=df['point_size'],
            sizemode='diameter',
            line=dict(width=0),
        ),
        hovertemplate='<br>'.join([
            'Chromosome: %{customdata[0]}',
            'Start: %{customdata[1]}',
            'Highlighted Gene(s): %{customdata[2]}',
            'Depth: %{customdata[3]}',
            'Weight: %{customdata[4]:.3f}',
            'Log2: %{y:.3f}',
        ]),
        customdata=df[['chromosome', 'start', 'gene', 'depth', 'weight']],
        name=name
    )

def create_gene_annotations(df: pd.DataFrame, color: str) -> List[Dict]:
    gene_groups = df.groupby('gene', group_keys=False).agg({
        'log2': lambda x: x.loc[x.abs().idxmax()],
        'genomic_position': 'first'
    }).reset_index()

    return [
        dict(
            x=row['genomic_position'],
            y=row['log2'],
            text=row['gene'],
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor=color,
            #ay=-60,
            font=dict(size=10, color=color),
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor=color,
            borderwidth=1,
            borderpad=4,
        )
        for _, row in gene_groups.iterrows()
    ]

def customize_layout(
    fig: go.Figure,
    df: pd.DataFrame,
    chrom_sizes: pd.DataFrame,
    chrom_centers: List[float],
    sample_name: str,
    annotations: Dict[str, List[Dict]],
    custom_genes: List[str],
    default_genes: List[str]
) -> None:
    """
    Customize the layout of the Plotly figure, including buttons for toggling gene annotations and tables.

    Args:
        fig: The Plotly figure to customize.
        df: The main DataFrame containing CNV data.
        chrom_sizes: DataFrame containing chromosome sizes.
        chrom_centers: List of chromosome center positions.
        sample_name: Name of the sample.
        annotations: Dictionary containing 'default' and 'custom' annotations.
        custom_genes: List of custom genes provided by the user.
        default_genes: List of default genes.
    """
    # Define Y-axis range with some padding
    y_min, y_max = df['log2'].min(), df['log2'].max()
    if abs(y_min) > y_max:
        y_padding = 0.1 * (abs(y_min))
        default_y_range = [y_min - y_padding, abs(y_min) + y_padding]
    else:
        y_padding = 0.1 * y_max
        default_y_range = [y_max - y_padding, y_max + y_padding]

    # Define Profile Y-axis range based on smoothed line
    profile_y_min, profile_y_max = df['smooth_log2'].min(), df['smooth_log2'].max()
    profile_padding = 0.1 * (profile_y_max - profile_y_min)
    profile_y_range = [profile_y_min - profile_padding, profile_y_max + profile_padding]

    # Define X-axis range
    x_min, x_max = df['genomic_position'].min(), df['genomic_position'].max()
    x_range = [x_min, x_max]

    # Define chromosome boundaries for vertical lines
    chrom_boundaries = chrom_sizes['cumulative_position'].tolist()

    # Combine annotations based on initial state (Show All Genes)
    initial_annotations = annotations['default'].copy()
    if custom_genes:
        initial_annotations += annotations['custom']

    # Prepare shapes for default and profile views
    shapes_default = []
    shapes_profile = []
    for boundary in chrom_boundaries[1:]:  # Skip the first boundary (start of chr1)
        # Shape for default view
        shapes_default.append(
            dict(
                type="line",
                x0=boundary,
                x1=boundary,
                y0=default_y_range[0],
                y1=default_y_range[1],
                line=dict(color="darkgrey", width=1.5, dash="solid"),
            )
        )
        # Shape for profile view
        shapes_profile.append(
            dict(
                type="line",
                x0=boundary,
                x1=boundary,
                y0=profile_y_range[0],
                y1=profile_y_range[1],
                line=dict(color="darkgrey", width=1.5, dash="solid"),
            )
        )

    # Initialize layout with default shapes
    fig.update_layout(
        title=dict(
            text=f"Genome-Wide Interactive CNV Plot - {sample_name}",
            font=dict(weight=1000, color='black'),
            xanchor='left',
            y=0.96,  # Move title higher
            yanchor='top'
        ),
        autosize=True,
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=50, r=50, t=50, b=50),
        hovermode='closest',
        annotations=initial_annotations,
        yaxis={'range': default_y_range},
        xaxis={'range': x_range},
        shapes=shapes_default  # Set initial shapes to default
    )

    # Update x-axis and y-axis for the main plot (first subplot)
    fig.update_xaxes(
        title="Chromosome",
        range=x_range,
        tickmode='array',
        tickvals=chrom_centers,
        ticktext=[chrom.replace('chr', '') for chrom in chrom_sizes['#chrom']],
        tickangle=0,
        tickfont=dict(size=10),
        showgrid=False,
        zeroline=True,
        zerolinecolor='black',
        zerolinewidth=1,
        row=1, col=1
    )

    fig.update_yaxes(
        title="Copy Number (Log2 ratio)",
        range=default_y_range,
        showgrid=True,
        gridcolor='lightgrey',
        gridwidth=0.05,
        zeroline=True,
        zerolinecolor='black',
        zerolinewidth=1,
        row=1, col=1
    )

    # Identify indices of highlighted traces and table traces
    highlighted_default_index = None
    highlighted_custom_index = None
    all_genes_table_index = None
    default_genes_table_index = None
    custom_genes_table_index = None
    empty_table_index = None

    for idx, trace in enumerate(fig.data):
        if trace.name == 'Highlighted Default Genes':
            highlighted_default_index = idx
        elif trace.name == 'Highlighted Custom Genes':
            highlighted_custom_index = idx
        elif trace.name == 'All Genes Table':
            all_genes_table_index = idx
        elif trace.name == 'Default Genes Table':
            default_genes_table_index = idx
        elif trace.name == 'Custom Genes Table':
            custom_genes_table_index = idx
        elif trace.name == 'Empty Table':
            empty_table_index = idx

    # Store table trace indices in a dictionary for easy access
    table_trace_indices = {
        "All Genes Table": all_genes_table_index,
        "Default Genes Table": default_genes_table_index,
        "Custom Genes Table": custom_genes_table_index,
        "Empty Table": empty_table_index
    }

    # Combine visibility flags for traces and tables
    # Initialize with main scatter and smoothed line visible
    visibility_base = [True, True]  # Traces 0 and 1

    # Append visibility for Highlighted Default and Custom Genes
    visibility_base.append(True if highlighted_default_index is not None else False)  # Trace 2: highlighted default genes
    visibility_base.append(True if highlighted_custom_index is not None else False)  # Trace 3: highlighted custom genes

    # Append visibility for tables
    visibility_base.append(True if all_genes_table_index is not None else False)    # Trace 4: All Genes Table
    visibility_base.append(False)  # Trace 5: Default Genes Table
    visibility_base.append(False)  # Trace 6: Custom Genes Table
    visibility_base.append(False)  # Trace 7: Empty Table

    # Now, set the initial visibility of all traces
    for i, trace in enumerate(fig.data):
        trace.visible = visibility_base[i]

    # Prepare annotations for different button states
    all_annotations = annotations['default'].copy()
    if custom_genes:
        all_annotations += annotations['custom']

    fig.update_layout(
        xaxis=dict(domain=[0, 1], anchor='y'),
        yaxis=dict(domain=[0.25, 0.9], anchor='x')  # Adjust this to move the plot down
    )

    # Adjust the position of the table
    fig.update_layout(
        xaxis2=dict(domain=[0, 1], anchor='y2'),
        yaxis2=dict(domain=[0, 0.2], anchor='x2')  # Adjust this to set table height
    )

    # Initialize buttons list with table trace indices
    buttons = create_buttons(
        fig,
        annotations,
        custom_genes,
        default_y_range,
        profile_y_range,
        table_trace_indices,
        shapes_default,
        shapes_profile
    )  # Pass fig here along with table_trace_indices and shapes

    # Update layout with updatemenus
    fig.update_layout(
        updatemenus=[{
            "buttons": buttons,
            "direction": "left",
            "pad": {"r": 10, "t": 10},
            "showactive": True,
            "x": 1.0,
            "xanchor": "right",
            "y": 1,
            "yanchor": "top"
        }],
    )

def initialize_visibility(traces):
    # Initialize all traces as invisible except the main scatter and smoothed line
    visibility = [False] * len(traces)
    visibility[0] = True  # Main scatter trace
    visibility[1] = True  # Smoothed line trace
    return visibility

def update_visibility(visibility: List[bool], trace_indices, state: bool) -> List[bool]:
    """
    Update the visibility of specified trace indices.

    Args:
        visibility (List[bool]): Current visibility list.
        trace_indices (int or List[int]): Single index or list of indices to update.
        state (bool): The visibility state to set.

    Returns:
        List[bool]: Updated visibility list.
    """
    if isinstance(trace_indices, int):
        visibility[trace_indices] = state
    elif isinstance(trace_indices, list):
        for idx in trace_indices:
            visibility[idx] = state
    return visibility



def create_button(label: str, visibility: List[bool], annotations: List[Dict], layout_updates: Dict = None) -> Dict:
    args = [{"visible": visibility}, {"annotations": annotations}]

    if layout_updates:
        args[1].update(layout_updates)

    return dict(
        args=args,
        label=label,
        method="update"
    )




def create_show_all_genes_button(
    visibility: List[bool],
    annotations: Dict[str, List[Dict]],
    default_y_range: List[float],
    table_trace_indices: Dict[str, int],
    shapes_default: List[Dict]
) -> Dict:
    new_visibility = visibility.copy()
    # Enable gene highlights
    new_visibility = update_visibility(new_visibility, [2, 3], True)  # Highlighted Default and Custom Genes
    # Set table visibility: Show All Genes Table, hide others
    new_visibility = update_visibility(new_visibility, [table_trace_indices["All Genes Table"]], True)
    new_visibility = update_visibility(new_visibility, [
        table_trace_indices["Default Genes Table"],
        table_trace_indices["Custom Genes Table"],
        table_trace_indices["Empty Table"]
    ], False)
    # Layout update for y-axis range and shapes
    layout_updates = {
        "yaxis.range": default_y_range,
        "shapes": shapes_default
    }
    return create_button(
        "Show All Genes",
        new_visibility,
        annotations['default'] + annotations['custom'],
        layout_updates
    )


def create_hide_all_genes_button(
    visibility: List[bool],
    default_y_range: List[float],
    table_trace_indices: Dict[str, int],
    shapes_default: List[Dict]
) -> Dict:
    new_visibility = visibility.copy()
    # Disable gene highlights
    new_visibility = update_visibility(new_visibility, [2, 3], False)  # Highlighted Default and Custom Genes
    # Disable smoothed line
    new_visibility = update_visibility(new_visibility, [1], False)
    # Set table visibility: Show Empty Table, hide others
    new_visibility = update_visibility(new_visibility, [table_trace_indices["Empty Table"]], True)
    new_visibility = update_visibility(new_visibility, [
        table_trace_indices["All Genes Table"],
        table_trace_indices["Default Genes Table"],
        table_trace_indices["Custom Genes Table"]
    ], False)
    # Layout update for y-axis range
    layout_updates = {
        "yaxis.range": default_y_range,
        "shapes": shapes_default
    }

    return create_button(
        "Hide All Genes",
        new_visibility,
        [],
        layout_updates
    )

def create_show_default_genes_button(
    visibility: List[bool],
    annotations: Dict[str, List[Dict]],
    default_y_range: List[float],
    table_trace_indices: Dict[str, int],
    shapes_default: List[Dict]
) -> Dict:
    new_visibility = visibility.copy()
    # Enable default genes, disable custom genes
    new_visibility = update_visibility(new_visibility, [2], True)  # Highlighted Default Genes
    new_visibility = update_visibility(new_visibility, [3], False)  # Highlighted Custom Genes
    # Set table visibility: Show Default Genes Table, hide others
    new_visibility = update_visibility(new_visibility, [table_trace_indices["Default Genes Table"]], True)
    new_visibility = update_visibility(new_visibility, [
        table_trace_indices["All Genes Table"],
        table_trace_indices["Custom Genes Table"],
        table_trace_indices["Empty Table"]
    ], False)
    # Layout update for y-axis range
    layout_updates = {
        "yaxis.range": default_y_range,
        "shapes": shapes_default
    }

    return create_button(
        "Show Default Genes",
        new_visibility,
        annotations['default'],
        layout_updates
    )


def create_show_custom_genes_button(
    visibility: List[bool],
    annotations: Dict[str, List[Dict]],
    default_y_range: List[float],
    table_trace_indices: Dict[str, int],
    shapes_default: List[Dict]
) -> Dict:
    new_visibility = visibility.copy()
    # Enable custom genes, disable default genes
    new_visibility = update_visibility(new_visibility, [3], True)  # Highlighted Custom Genes
    new_visibility = update_visibility(new_visibility, [2], False)  # Highlighted Default Genes
    # Set table visibility: Show Custom Genes Table, hide others
    new_visibility = update_visibility(new_visibility, [table_trace_indices["Custom Genes Table"]], True)
    new_visibility = update_visibility(new_visibility, [
        table_trace_indices["All Genes Table"],
        table_trace_indices["Default Genes Table"],
        table_trace_indices["Empty Table"]
    ], False)
    # Layout update for y-axis range
    layout_updates = {
        "yaxis.range": default_y_range,
        "shapes": shapes_default
    }

    return create_button(
        "Show Custom Genes",
        new_visibility,
        annotations['custom'],
        layout_updates
    )


def create_show_profile_button(
    visibility: List[bool],
    default_y_range: List[float],
    profile_y_range: List[float],
    existing_annotations: List[Dict],
    table_trace_indices: Dict[str, int],
    shapes_profile: List[Dict]
) -> Dict:
    new_visibility = visibility.copy()
    # Disable main scatter and gene highlights
    new_visibility = update_visibility(new_visibility, [0, 2, 3], False)
    # Enable smoothed line
    new_visibility[1] = True
    # Set table visibility: Show Empty Table, hide others
    new_visibility = update_visibility(new_visibility, [table_trace_indices["Empty Table"]], True)
    new_visibility = update_visibility(new_visibility, [
        table_trace_indices["All Genes Table"],
        table_trace_indices["Default Genes Table"],
        table_trace_indices["Custom Genes Table"]
    ], False)
    # Layout update for y-axis range
    layout_updates = {
        "yaxis.range": profile_y_range,
        "shapes": shapes_profile
    }
    # Clear annotations when showing profile
    current_annotations = []  # Assuming you want to hide annotations when showing profile
    return create_button(
        "Show Profile",
        new_visibility,
        current_annotations,
        layout_updates
    )



def create_buttons(
    fig: go.Figure,
    annotations: Dict[str, List[Dict]],
    custom_genes: List[str],
    default_y_range: List[float],
    profile_y_range: List[float],
    table_trace_indices: Dict[str, int],
    shapes_default: List[Dict],
    shapes_profile: List[Dict]
) -> List[Dict]:
    visibility = initialize_visibility(fig.data)
    buttons = [
        create_show_all_genes_button(
            visibility,
            annotations,
            default_y_range,
            table_trace_indices,
            shapes_default
        ),
        create_show_default_genes_button(
            visibility,
            annotations,
            default_y_range,
            table_trace_indices,
            shapes_default
        ),
    ]

    if custom_genes:
        buttons.append(
            create_show_custom_genes_button(
                visibility,
                annotations,
                default_y_range,
                table_trace_indices,
                shapes_default
            )
        )

    buttons.append(
        create_hide_all_genes_button(
            visibility,
            default_y_range,
            table_trace_indices,
            shapes_default
        )
    )

    buttons.append(
        create_show_profile_button(
            visibility,
            default_y_range,
            profile_y_range,
            annotations['default'] + annotations['custom'],
            table_trace_indices,
            shapes_profile
        )
    )

    return buttons



def create_gene_table(df: pd.DataFrame, genes: List[str], name: str) -> go.Table:
    """
    Creates a table with data for selected genes.
    """
    gene_data = df[df['gene'].isin(genes)]  # Removed .copy() as it's unnecessary
    gene_data = gene_data.sort_values(['gene', 'start'])

    table = go.Table(
        header=dict(
            values=['Gene', 'Chromosome', 'Start', 'Log2', 'Depth', 'Weight'],
            font=dict(size=12, weight=1000, color='white'),
            align='left',
            height=30,
            fill_color='darkblue',
        ),
        cells=dict(
            values=[
                gene_data['gene'],
                gene_data['chromosome'],
                gene_data['start'],
                gene_data['log2'].round(3),
                gene_data['depth'],
                gene_data['weight'].round(3)
            ],
            font=dict(size=11),
            align='left',
            height=25
        ),
        columnwidth=[1, 1, 1, 1, 1, 1],  # Equal width for all columns
        name=name,
        visible=False  # Initially set to invisible
    )

    return table

def process_sample(
    cnv_file: str,
    chrom_sizes: pd.DataFrame,
    output_dir: str,
    prefix: str,
    args_genes: str
):
    try:
        df = read_data(cnv_file)
        if df.empty:
            logging.warning(f"Skipping empty dataset: {cnv_file}")
            return

        print(f"Processing {cnv_file}")

        df = calculate_genomic_position(df, chrom_sizes)
        df = df[df['chromosome'].isin(CLASSIC_CHROMOSOMES)]

        default_genes, custom_genes = get_highlight_genes(args_genes)
        all_genes = default_genes + custom_genes  # Concatenate once and reuse

        # Create scatter plot traces and annotations
        scatter_traces, annotations = create_scatter_plot(df, default_genes, custom_genes)

        # Create figure with subplots: main plot and table
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.8, 0.2],
            specs=[[{"type": "scatter"}],
                   [{"type": "table"}]],
            vertical_spacing=0.1
        )

        # Add scatter traces to the main plot
        for trace in scatter_traces:
            fig.add_trace(trace, row=1, col=1)

        # Add gene tables to the subplot
        fig.add_trace(create_gene_table(df, all_genes, "All Genes Table"), row=2, col=1)
        fig.add_trace(create_gene_table(df, default_genes, "Default Genes Table"), row=2, col=1)
        fig.add_trace(create_gene_table(df, custom_genes, "Custom Genes Table"), row=2, col=1)
        fig.add_trace(create_gene_table(df, [], "Empty Table"), row=2, col=1)

        sample_name = os.path.basename(cnv_file).split('.')[0]
        customize_layout(
            fig,
            df,
            chrom_sizes,
            chrom_sizes['chrom_center'].tolist(),
            sample_name,
            annotations,
            custom_genes,
            default_genes
        )

        output_file = os.path.join(output_dir, f"{prefix}_{sample_name}.html")
        config = {
            'modeBarButtonsToAdd': ['hoverclosest'],
            'displayModeBar': True,
            'displaylogo': True,
            'scrollZoom': True,
            'toImageButtonOptions': {
                'format': 'png', # one of png, svg, jpeg, webp
                'filename': f'image_{sample_name}',
                'height': 900,
                'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
            }
        }

        fig.write_html(
            output_file,
            include_plotlyjs='cdn',
            config=config
        )

        print(f"Plot saved to {output_file}")
    except Exception as e:
        logging.error(f"Failed to process {cnv_file}: {e}")
        # Optionally, continue with the next file or task

def main():
    setup_logging()
    args = parse_arguments()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    logging.info("Loading chromosome sizes data")
    chrom_sizes = read_data(args.chrom_sizes, is_chrom_sizes=True)

    for filename in os.listdir(args.input_dir):
        if filename.endswith('.cnr'):
            cnv_file = os.path.join(args.input_dir, filename)
            logging.info(f"Processing CNV file: {cnv_file}")
            process_sample(cnv_file, chrom_sizes, args.output_dir, args.prefix, args.genes)

    logging.info("All samples processed.")
    logging.info("Application completed successfully.")

if __name__ == "__main__":
    main()