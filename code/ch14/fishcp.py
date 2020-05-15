from bokeh.plotting import figure, show, output_file
from bokeh.models.markers import marker_types
from bokeh.transform import factor_cmap, factor_mark
from pandas import read_csv

df = read_csv('../samples/fishdata.csv')
df['feeds_and_species'] = df['feeds'] + ', ' + df['species']
all_markers = [mt for mt in marker_types]
SPECIES = list(set(df['species']))
MARKERS = all_markers[:len(SPECIES)]
feeds = list(set(df['feeds']))
ttl = 'Metabolic variations based on 1H NMR profiling of fishes'

p = figure(plot_height=600, plot_width=700, title = ttl)
p.xaxis.axis_label = 'Principal Component 1: 35.8%'
p.yaxis.axis_label = 'Principal Component 2: 15.1%'
p.scatter('PC1', 'PC2', source=df, size=12, fill_alpha=0.3, 
          marker=factor_mark('species', MARKERS, SPECIES),
          color=factor_cmap('feeds', 'Category10_3', feeds),
          legend_field='feeds_and_species')
p.legend.location = 'top_left'
p.legend.click_policy = 'hide'
output_file('scatter.html')
show(p)
