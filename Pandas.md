Set a column based on a pd.query result:
>df.loc[df.query('ColumnA > 3').index, "new_column"] = "Triplet"
