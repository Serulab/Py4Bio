from bokeh.charts import Scatter, output_file, show

x = [1, 2, 3, 4, 5, 6, 7, 8]
y = [2.1, 6.45, 3, 1.4, 4.55, 3.85, 5.2, 0.7]
z = [.5, 1.1, 1.9, 2.5, 3.1, 3.9, 4.85, 5.2]
species = ['cat', 'cat', 'cat', 'dog', 'dog', 'dog', 'mouse', 'mouse']
country = ['US', 'US', 'US', 'US', 'UK', 'UK', 'BR', 'BR']

df = {'time': x, 'weight 1': y, 'weight 2': z, 'species':species, 'country': country}

scatter = Scatter(df, x='time', y='weight 1', color='country', marker='species',
                  title="Auto MPG", xlabel="Time in days",
                  ylabel="Weight in grams")

output_file('scatter.html')
show(scatter)
