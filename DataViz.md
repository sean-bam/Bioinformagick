layout: page
title: "Data Visualization with Plotly"
permalink: /DataViz

## <a name='TOC'>Table of Contents</a>

1. [Sort rows of a dataframe by UPGMA](#UPGMA)
2. [Build a heatmap](#Heatmap)
3. [Scatter plots with shared axis](#ScatterShared)


**<a name="UPGMA">Sort rows of a dataframe by UPGMA.</a>**

Plotly has a thin-wrapper around SciPy's hierarchical clustering functions. The default is to use UPGMA. Using a dataframe of distances, we can run UPGMA clustering through plotly, and use the resulting leaf order to sort the dataframe.
```python
fig = ff.create_dendrogram(df, labels = df.index.tolist())
dendro_leaves = fig['layout']['xaxis']['ticktext']
df2 = df.loc[dendro_leaves,:]
```

**<a name="Heatmap">Build a heatmap.</a>**
```python
fig = go.Figure(data=go.Heatmap(
                    x=df.columns,
                    y=df.index,
                    z=df.values,
                    colorscale = 'Turbo',
                    #colorscale = [
                    #    # Colors for values of zero
                    #    [0, "rgb(0,0,0)"],
                    #    # Colors for values b/n 0.1-20% of max 
                    #    [0.2, "rgb(0,0,0)"],
                    #    # Colors for values b/n 20-40% of max
                    #    [0.4, "rgb(0,250,253)"],
                    #    # Colors for values b/n 40-60% of max
                    #    [0.6, "rgb(3,255,11)"],
                    #    # Colors for values b/n 60-80% of max
                     #   [0.8, "rgb(245,255,11)"],
                     #   # Colors for values b/n 80-100% of max
                     #   [1, "rgb(243,20,21)"]],
                    #zmin = 0.00001,
                    zmax = 1,
                    #zsmooth = 'best',
                    colorbar_thickness = 15,
                    #colorbar_tick0 = 0,
                    #colorbar_dtick = 1,
                    colorbar_ticks = "outside"))
    
#update X and Y axes
fig.update_yaxes(title_text = 'Y-axis-title',
                    titlefont=dict(size=16,
                                  family='Arial',
                                  color='black'
                                  ),
                    showgrid=False,
                    automargin = True,
                    showticklabels = True,
                    showline=True,
                    mirror=True,
                    linewidth=0.5, 
                    linecolor='black',
                )

fig.update_xaxes(title_text = 'X-axis-title',
                    titlefont=dict(size=16,
                                  family='Arial',
                                  color='black'
                                  ),
                    showgrid=False,
                    automargin = True,
                    showticklabels = False,
                    showline=True,
                    #mirror=True,
                    linewidth=0.5, 
                    linecolor='black',
                )

fig.update_layout(title_text = 'Figure title',
                    titlefont = dict(size=16,
                                    family='Arial', 
                                    color='Black'
                                    )
                  )
#fig.show()
#fig.write_html('heatmap.html', default_height = 800)
#fig.write_image('path/to/svg',  height = 1200, width = 800)
````

**<a name="ScatterShared">Scatter plots with shared axis</a>**
```python
fig = make_subplots(rows=3, 
                    cols=1,
                    shared_xaxes = True
                   )

fig.add_trace(go.Scatter(x=df.x.tolist(),
                         y=df.y.tolist()
                        ), 
              row=1, 
              col=1)

fig.add_trace(go.Scatter(x=df2.x.tolist(), 
                         y=df2.y.tolist()
                        ), 
              row=2, 
              col=1)

fig.add_trace(go.Scatter(x=df3.x.tolist(), 
                         y=df3.y.tolist()
                        ), 
              row=3, 
              col=1)


#this doesn't show in jupyter for some reason, but does show up in the SVG file
fig.add_vline(x=140,
             line_dash="dot")

#update axes of plots
fig.update_yaxes(title_text="ytitle1", row=1, col=1)
fig.update_yaxes(title_text="ytitle2", row=2, col=1)
fig.update_yaxes(title_text="ytitle3", row=3, col=1)
fig.update_xaxes(title_text="xtitle", row=4, col=1)

#update figure size
fig.update_layout(width=600, 
                  height=600, 
                  template= 'simple_white'
                  )
#fig.show()
#fig.write_image('path/to/svg')
```
