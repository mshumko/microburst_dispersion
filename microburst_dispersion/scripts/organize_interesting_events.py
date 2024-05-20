"""
Organize interesting microburst events into folders separated by category type.
"""
import pathlib

import pandas as pd

import microburst_dispersion

interesting_microbust_list = pd.read_excel(
    'C:\\Users\\shumkms1\\Documents\\research\\microburst_dispersion\\data\\'
    'firebird_interesting_dispersion_examples.xlsx',
    skiprows=1, 
    usecols='A:C'  # We don't need any of the other columns for now.
    )

input_plot_dir = save_dir = microburst_dispersion.config['here'].parent / 'plots' / 'dispersion_summary'
output_plot_dir = input_plot_dir / 'interesting'

if not pathlib.Path(output_plot_dir).exists():
    output_plot_dir.mkdir()

for i, row in interesting_microbust_list.iterrows():
    input_plot_path = input_plot_dir / (row['plot_filename'] + '.png')

    # This is how you copy using pathlib.
    output_name = (row['plot_filename'] + '.png')
    output_name_chunks = output_name.split('_')
    output_name_chunks[3] = f'{row["dispersion"]}_{output_name_chunks[3]}'
    output_plot_path = output_plot_dir / '_'.join(output_name_chunks)
    output_plot_path.write_bytes(input_plot_path.read_bytes())

pass