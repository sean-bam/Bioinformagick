def parse_dali(dalioutput):
    """
    Gets the top 9 hits from a dali 'summary' output
    Clearly, needs improvement..
    """
    data = []
    with open(dalioutput) as f:
        for line in f:
            if line.startswith("# Query:"):
                query = line.strip().split()[2]
            if line.startswith("   "):
                hit_no, chain, z, rmsd, lali, nres, identity, description = line.strip().split(maxsplit=7)
                data.append({"query_id" : query,
                             "chain" : chain,
                             "z" : z,
                             "rmsd" : rmsd,
                             "lali" : lali,
                             "nres" : nres,
                             "identity" : identity,
                             "description" : description})
    df = pd.DataFrame.from_records(data)
    return df