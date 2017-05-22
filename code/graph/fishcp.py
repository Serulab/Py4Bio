from bokeh.charts import Scatter, output_file, show
from pandas import DataFrame

df = DataFrame.from_csv('fishdata.csv')

scatter = Scatter(df, x='PC1', y='PC2', color='feeds',
        marker='species', title=
        'Metabolic variations based on 1H NMR profiling of fishes',
        xlabel='Principal Component 1: 35.8%',
        ylabel='Principal Component 2: 15.1%')
scatter.legend.background_fill_alpha = 0.3
output_file('scatter.html')
show(scatter)
