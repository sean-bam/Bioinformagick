Sort rows by UPGMA
```
fig = ff.create_dendrogram(df_pathway2, labels = df_pathway2.index.tolist())
dendro_leaves = fig['layout']['xaxis']['ticktext']
df_pathway3 = df_pathway2.loc[dendro_leaves,:]
```

Make a heatmap
```
fig = go.Figure(data=go.Heatmap(
                    x=df_pathway4.columns,
                    y=df_pathway4.index,
                    z=df_pathway4.values,
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
fig.update_yaxes(title_text = 'Pathway',
                    titlefont=dict(size=16,family='Arial', color='black'),
                    showgrid=False,
                    automargin = True,
                    showticklabels = True,
                    showline=True,
                    mirror=True,
                    linewidth=0.5, 
                    linecolor='black')

fig.update_xaxes(title_text = 'Transposons',
                    titlefont=dict(size=16,family='Arial', color='black'),
                    showgrid=False,
                    automargin = True,
                    showticklabels = False,
                    showline=True,
                    #mirror=True,
                    linewidth=0.5, 
                    linecolor='black')

fig.update_layout(title_text = 'Transposon-mobilized COG pathways',
                    titlefont = dict(size=16,family='Arial', color='Black'))
fig.show()
#fig.write_html('COG_pathway_heatmap.html', default_height = 800)
#fig.write_image('../99figures/COG_pathway_heatmap2.svg',  height = 1200, width = 800)
````

Scatter plots with shared xaxis
```
fig = make_subplots(rows=3, 
                    cols=1,
                    shared_xaxes = True
                   )

fig.add_trace(go.Scatter(x=df_params.resolution.tolist(),
                         y=df_params.no_subcommunities.tolist()
                        ), 
              row=1, 
              col=1)

fig.add_trace(go.Scatter(x=df_params.resolution.tolist(), 
                         y=df_params.community_size.tolist()
                        ), 
              row=2, 
              col=1)

fig.add_trace(go.Scatter(x=df_params.resolution.tolist(), 
                         y=df_params.persistence.tolist()
                        ), 
              row=3, 
              col=1)


#this doesn't show in jupyter for some reason, but does show up in the SVG file
fig.add_vline(x=140,
             line_dash="dot")

#update axes of plots
fig.update_yaxes(title_text="# communities", row=1, col=1)
fig.update_yaxes(title_text="Community size", row=2, col=1)
fig.update_yaxes(title_text="Persistence", row=3, col=1)
fig.update_xaxes(title_text="Resolution", row=4, col=1)

#update figure size
fig.update_layout(width=600, height=600, template= 'simple_white')
#fig.show()
fig.write_image('hidef_params.svg')
```
