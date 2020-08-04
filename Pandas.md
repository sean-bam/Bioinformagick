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
