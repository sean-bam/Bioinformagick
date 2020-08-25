Set a column based on a pd.query result:
>df.loc[df.query('ColumnA > 3').index, "new_column"] = "Triplet"

Explode a cell into rows ([Source](https://stackoverflow.com/questions/17116814/pandas-how-do-i-split-text-in-a-column-into-multiple-rows/21032532))
```
s = df['Members'].str.split(',').apply(pd.Series, 1).stack()   #Split the cells in "Members" by commas
s.index = s.index.droplevel(-1)                                # to line up with df's index
s.name = 'Members'                                             # needs a name to join
del df['Members']
df2 = df.join(s)
```

Set a column with the count of elements in another column.
```
df["size"] = (df.groupby(by = ['mash_cluster_rep'])['mash_cluster_rep']
                                .transform('count')
                             )
```

Update values in columns of one dataframe using columns of another dataframe. Only matching columns are updated!
[Source](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.update.html)
```
df1.update(df)
```
