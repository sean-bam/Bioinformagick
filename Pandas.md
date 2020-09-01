Set a new column based on a `query` result:
>df.loc[df.query('ColumnA > 3').index, "new_column"] = "Triplet"

Update a column based on multiple conditions using `mask` or `where`. Note: Use `&` instead of `and`.
```
#Make boolean criteria
cond1 = df.Genus == "-"
cond2 = df.Family == "Siphoviridae"
cond3 = df.Family == "Myoviridae"
cond4 = df.Family == "Podoviridae"
cond5 = df.Subfamily.isna()

#Mask the dataframe using the boolean criteria, which preserves the shape of the dataframe. 
#Unmasked values become NaN, so fill those
#return the new column
df["Family"] = (df.mask(cond1 & (cond2 | cond3 | cond4) & cond5)
                                .fillna({"Family" : "-"})
                                .loc[:,"Family"])

```

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
#Option1
df["size"] = (df.groupby(by = ['mash_cluster_rep'])['mash_cluster_rep']
                                .transform('count')
                             )
                             
#Option2
df = (pd.read_csv('old/pevzner_candidate_phage_annotated.csv')
        .fillna('None') #groupby drops the row if NaN (?)
        .groupby(['Order','Family','Subfamily', 'Genus', "NodeID"])
        .size()
        .reset_index(name='Count')
      )
```

Update values in columns of one dataframe using columns of another dataframe. Only matching columns are updated!
[Source](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.update.html)
```
df1.update(df)
```
