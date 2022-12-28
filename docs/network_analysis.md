
Build a networkx object from a dataframe, check density
```
import networkx as nx
G=nx.from_pandas_adjacency(df_adj)

#check density
nx.density(G)
```

