**Set a new column based on a `query` result:**
```
df.loc[df.query('ColumnA > 3').index, "new_column"] = "Triplet"

WARNING: This produces some unexpected behavior if you do something like:
list_of_identifiers = ['abc', 'def']
df.loc[df.query('ColumnA in @list_of_identifiers').index, "new_column"] = "Triplet"

Sometimes, values that ARENT in list_of_identifiers are set

So, if using a list, use mask/where approach described below
```


**update a column using `where` or 'mask'**
```
#Change all values in the 'icity' column that are greater than 1 to 1
df["icity"] = df.where(df.icity <1, 1).icity

#the same, but using mask
df.mask(df.icity >1)
df["icity"] = df.mask(df.icity >1, 1).icity

#the same, but using a list
list_of_identifiers = ['abc', 'def']
df["icity"] = np.nan
df["icity"] = df.mask(df.icity.isin(list_of_identifiers), "new_value").icity
```

**Update a column based on multiple conditions using `mask`**
Note: Use `&` instead of `and`.
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
**Update 2+ columns based on a conditions using `mask`**
For example, change a dataframe where some rows contain start > stop
to start > stop
```
df[["start","stop"]] = (df.loc[:,["start","stop"]]          #select the columns to update
                          .mask(df.start > df.stop,         #select rows where condition is true
                                df[["stop","start"]],       #replace by swapping columns
                                axis = 'index')
                             )
                         )
```


**Explode a cell into rows** ([Source](https://stackoverflow.com/questions/17116814/pandas-how-do-i-split-text-in-a-column-into-multiple-rows/21032532))
```
s = df['Members'].str.split(',').apply(pd.Series, 1).stack()   #Split the cells in "Members" by commas
s.index = s.index.droplevel(-1)                                # to line up with df's index
s.name = 'Members'                                             # needs a name to join
del df['Members']
df2 = df.join(s)
```

**Set a column with the count of elements in another column.**
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
**How to Split-Apply-Combine**
```
#Define a function you want to apply to the group
def calculate_length(group):
    group["length"] = group["end"] - group["start"]
    
#Split the dataframe by something, apply the function. The results are automatically applied to every group
df2 = df.groupby('accession').apply(calculate_length)

#Rest the index to return the dataframe to its original shape
df3 = df2.reset_index(drop = True)
```

**Update values in columns of one dataframe using columns of another dataframe.** 
Note: Only matching indexes+columns are updated! Use `df.set_index()` to change indexes
Note: Duplicate indexes mess this up
[Source](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.update.html)
```
df1.update(df)
```

**Select rows before and after a given row**
[Source](https://stackoverflow.com/questions/48630060/select-n-rows-above-and-below-a-specific-row-in-pandas)
```
#Select rows of interest
idxs = df.query('name == "Mn_catalase_like" or name == "Dps"').index

#Loop through the list, slice the dataframe using iloc and the target index, +/- 2 rows
df_list = []
for idx in idxs:
    df1 = df.iloc[idx - 2 : idx + 2]
    df_list.append(df1)
    
#combine the results into a single dataframe
df2 = pd.concat(df_list)
```

**Split a dataframe into groups, set a column with the results of a named aggregation**
```
grouped = df_psiblast4.groupby(['phage', 'contig'])

#Make a column named "shared_orfs" that takes the column "protein_id" and calculates its size.
df_shared = grouped.agg(shared_orfs=pd.NamedAgg(column='protein_id', 
                                                   aggfunc='size'
                                                  )
                          )

df_shared2 = df_shared.reset_index()
```

**Remove reciprocal entries**
```
df = df.loc[:,["qaccver", "q_ign", "subj_ign", "dist", "IR_dist"]]

#create a column that joins Seq1 - Seq2 or Seq2 - Seq1 to Seq1Seq2
df["pairs"] = df.apply(lambda row: ''.join(sorted([row["q_ign"], row["subj_ign"]])), axis = 1)

#remove rows with no matching pair and sort the database
df2 = df.drop_duplicates(subset = 'pairs').drop(columns = ['pairs'])
```

**Get the first group in a groupby object**
```
group = df.groupby('motif_name')
df2 = group.get_group((list(group.groups)[0]))
```

**Bin continuous data into custom bins**
From [here](https://towardsdatascience.com/histograms-with-plotly-express-complete-guide-d483656c5ad7)
```
bins = [0,25,50,75,100,1000,10000,100000]
df['counts'] = pd.cut(df['dist_from_att'], bins=bins, include_lowest = True)
df_hist = df["counts"].value_counts().sort_index().to_frame().reset_index()
df_hist.rename(columns = {'index' : 'bins'}, inplace = True)
df_hist["bins"] = df_hist["bins"].astype('str')
```

**Get the most common value of a group**
The trick is that you have to use a lambda function to get the mode, because there may be multiple values, so this function just gets the first
```
df["top_cluster"] = df.groupby('tnsb_node')['cluster'].transform(lambda x: x.mode()[0])
```

**Remove duplicate columns based on the columns 'protein' and 'profile', but ignoring Na's in 'protein'**
[link](https://stackoverflow.com/questions/50154835/drop-duplicates-but-ignore-nulls)
```
df = df[(~df[['protein', 'profile']].duplicated()) | df['protein'].isna()]
```

**Create a dataframe from a for loop**
Iteratively adding things to a dataframe is slower than creating a list/dictionary, then creating a dataframe from that. 
Here is one example:
```
data = []
    with open(islandfile) as f:
        for line in f:
            if line.startswith("==="):
                drop1, accession, num_features, assembly, accession2, other = line.strip().split(maxsplit = 5)
                data.append({"accession" : accession,
                             "num_features" : num_features,
                             "assembly" : assembly,
                             "accession2" : accession2,
                             "other" : other})
df = pd.DataFrame.from_records(data)
```

**Select rows where there is at least one element in the specified column(s)**
```
df.loc[df.loc[:, ["protein_id", "locus_tag", "product"]].any(axis='columns')]
```

**Convert two columns from a dataframe into a python dictionary**
```
my_dict = dict(df.loc[:,["keycolumn","valuecolumn"]].values)
```

**remove spaces in column names**
```
columns = []
for col in df.columns.tolist():
    columns.append(col[1].replace(" ", "_"))
df.columns = columns
```

