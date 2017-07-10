from bokeh.charts import HeatMap, bins, output_file, show
import pandas as pd
DATA_FILE = '../../samples/GSM188012.CEL'
dtype = {'x': int, 'y': int, 'lux': float}
dataset = pd.read_csv(DATA_FILE, sep='\t', dtype=dtype)
hm = HeatMap(dataset, x=bins('x'), y=bins('y'), values='lux',
            title='Expression', stat='mean')
output_file("heatmap7.html", title="heatmap.py example")
show(hm)
