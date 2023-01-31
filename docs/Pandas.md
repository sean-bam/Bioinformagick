layout: page
title: "Pandas"
permalink: /Pandas

## <a name='TOC'>Pandas</a>

1. [Set a new column based on a `query` result](#newcolumnusingquery)
2. [Update a column using `where` or `mask`](#updatecolumnusingwhere)
3. [Update a column based on multiple conditions using `mask`](#Updateacolumnusingmultiple)
4. [Update multiple columns using `mask`](#Update2cols)
5. [Explode a cell into rows](#explode)
6. [Set a column with the count of elements in another column](#Countelementswithgroupby)
7. [How to split-apply-combine](#SplitApplyCombine)
8. [Update the values of one dataframe using another](#UseUpdate)
9. [Select rows before/after a row of interest](#SelectBeforeAfter)
10. [Set a new column using named aggregation](#UseNamedAggregation)
11. [Remove reciprocal entries/duplicates](#RemoveReciprocal)
12. [Get the first group in a `groupby` dataframe](#GetFirstGroup)
13. [Bin continuous data and visualize with a histogram](#BinContinuous)
14. [Get the mode of a `groupby` dataframe](#GetGroupMode)
15. [Remove duplicates but ignore Na's](#RemoveDupsIgnoringNA)
16. [Construct a dataframe using a `for` loop](#CreateDataframeFromForLoop)
17. [Select rows with at least one non-Na in the given column(s)](#SelectRowsWithOneElement)
18. [Convert a dataframe into a dictionary](#DataframeToDict)
19. [Remove spaces in column names](#RemoveSpacesinColNames)

**<a name="newcolumnusingquery">Set a new column based on a `query` result:</a>**

```python
df.loc[df.query('ColumnA > 3').index, "new_column"] = "Triplet"
```

**<a name="updatecolumnusingwhere">update a column using `where` or `mask`:</a>**
```python
import numpy as np
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

**<a name="Updateacolumnusingmultiple">Update a column based on multiple conditions using `mask`</a>**

Note: Use `&` instead of `and`.
```python
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
<a name="Update2cols">**Update 2+ columns based on a conditions using `mask`</a>**

For example, change a dataframe where some rows contain start > stop to start < stop
```python
df[["start","stop"]] = (df.loc[:,["start","stop"]]          #select the columns to update
                          .mask(df.start > df.stop,         #select rows where condition is true
                                df[["stop","start"]],       #replace by swapping columns
                                axis = 'index')
                             )
                         )
```


<a name="explode">**Explode a cell into rows**</a>
```python
df["profile"] = df.profile.str.split(",").tolist()
df = df.explode('profile')
```

<a name="Countelementswithgroupby">**Set a column with the count of elements in another column.**</a>

Note: `groupby` drops the row if NaN. So, a second option is shown to deal with this
```python
#Option1
df["size"] = (df.groupby(by = ['mash_cluster_rep'])['mash_cluster_rep']
                                .transform('count')
                             )
                             
#Option2
df = (pd.read_csv('path/to/csv')
        .fillna('None') 
        .groupby(['Family'])
        .size()
        .reset_index(name='Count')
      )
```

<a name="SplitApplyCombine">**How to Split-Apply-Combine**</a>
```python
#Define a function you want to apply to the group
def calculate_length(group):
    group["length"] = group["end"] - group["start"]
    return group
    
#Split the dataframe by something, apply the function. The results are automatically applied to every group
df2 = df.groupby('accession').apply(calculate_length)

#Rest the index to return the dataframe to its original shape
df3 = df2.reset_index(drop = True)
```

<a name="UseUpdate">**Update values in columns of one dataframe using columns of another dataframe.**</a>

Note: Only matching indexes+columns are updated! Use `df.set_index()` to change indexes

Note: Duplicate indexes mess this up
[Source](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.update.html)
```python
df2.update(df)
```

<a name="SelectBeforeAfter">**Select rows before and after a given row**</a>
[Source](https://stackoverflow.com/questions/48630060/select-n-rows-above-and-below-a-specific-row-in-pandas)
```python
#Select rows of interest
idxs = df.query('name == "something" or name == "something_else"').index

#Loop through the list, slice the dataframe using iloc and the target index, +/- 2 rows
df_list = []
for idx in idxs:
    df1 = df.iloc[idx - 2 : idx + 2]
    df_list.append(df1)
    
#combine the results into a single dataframe
df2 = pd.concat(df_list)
```
Option 2:
```python
#assuming the dataframe has a column "seed" that is set to True
#get 3 rows around these rows

def get_minimum_from_seed(group):
  all_indexes = group.index
  seed_indexes = group.query('seed = True').index
  group["dist_from_seed"] = [minimum(index, seed_indexes) for index in all_indexes)
  
  return group

df2 = df.groupby('contig').apply(get_minimum_from_seed)
df3 = df2.query('dist_from_seed < 3')
```

<a name="UseNamedAggregation">**Split a dataframe into groups, set a column with the results of a named aggregation**</a>
```python
grouped = df.groupby(['contig'])

#Make a column named "shared_orfs" that takes the column "protein_id" and calculates its size.
df_shared = grouped.agg(shared_orfs=pd.NamedAgg(column='protein_id', 
                                                   aggfunc='size'
                                                  )
                          )

df_shared2 = df_shared.reset_index()
```

<a name="RemoveReciprocal">**Remove reciprocal entries**</a>
```python
#create a column that joins Query : Subject or Subject : Query to QuerySubject
df["pairs"] = df.apply(lambda row: ''.join(sorted([row["query"], row["subject"]])), axis = 1)

#remove rows with no matching pair and sort the database
df2 = df.drop_duplicates(subset = 'pairs').drop(columns = ['pairs'])
```

<a name="GetFirstGroup">**Get the first group in a groupby object**</a>
```python
group = df.groupby('motif_name')
df2 = group.get_group((list(group.groups)[0]))
```

<a name="BinContinuous">**Bin continuous data into custom bins**</a>
From [here](https://towardsdatascience.com/histograms-with-plotly-express-complete-guide-d483656c5ad7).
Also shown is a way to visualize the results as a histogram
```python
bins = [0,25,50,75,100,1000,10000,100000]
df['counts'] = pd.cut(
                    df['dist_from_att'], 
                    bins=bins, 
                    include_lowest = True
                    )
df_hist = (df["counts"]
            .value_counts()
            .sort_index()
            .to_frame()
            .reset_index()
            )
            
df_hist.rename(columns = {'index' : 'bins'}, inplace = True)
df_hist["bins"] = df_hist["bins"].astype('str')
```

<a name="GetGroupMode">**Get the mode from a group**</a>
The trick is that you have to use a lambda function to get the mode, because there may be multiple values, so this function just gets the first
```python
df["top_cluster"] = df.groupby('tnsb_node')['cluster'].transform(lambda x: x.mode()[0])
```

<a name="RemoveDupsIgnoringNA">**Remove duplicate columns based on the columns 'protein' and 'profile', but ignoring Na's in 'protein'**</a>
[link](https://stackoverflow.com/questions/50154835/drop-duplicates-but-ignore-nulls)
```python
df = df[(~df[['protein', 'profile']].duplicated()) | df['protein'].isna()]
```

<a name="CreateDataframeFromForLoop">**Create a dataframe from a for loop**</a>
Iteratively adding things to a dataframe is slower than creating a list/dictionary, then creating a dataframe from that. 
Here is one example:
```python
data = []
    with open(islandfile) as f:
        for line in f:
            if line.startswith("==="):
                drop1, accession, num_features, assembly, accession2, other = line.strip().split(maxsplit = 5)
                data.append({"accession" : accession,
                             "num_features" : num_features,
                             "assembly" : assembly,
                             "accession2" : accession2,
                             "other" : other,
                             }
                           )
df = pd.DataFrame.from_records(data)
```

<a name="SelectRowsWithOneElement">**Select rows where there is at least one element in the specified column(s)**</a>
```python
df.loc[df.loc[:, ["protein_id", "locus_tag", "product"]].any(axis='columns')]
```

<a name="DataframeToDict">**Convert two columns from a dataframe into a python dictionary**</a>
```python
my_dict = dict(df.loc[:,["keycolumn","valuecolumn"]].values)
```

<a name="RemoveSpacesinColNames">**remove spaces in column names**</a>
```python
columns = []
for col in df.columns.tolist():
    columns.append(col[1].replace(" ", "_"))
df.columns = columns
```

