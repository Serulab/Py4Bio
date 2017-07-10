from bokeh.charts import output_file, Chord
from bokeh.io import show
import pandas as pd
data = pd.read_csv('../../samples/test3.csv')
chord_from_df = Chord(data, source='name_x', target='name_y',
                      value='value')
output_file('chord.html')
show(chord_from_df)
