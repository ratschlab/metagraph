import pandas as pd


def df_from_search_result(json_res):
    def _build_dict(row, result):
        d = dict(row)
        if 'properties' in d.keys():
            props = d.pop('properties')
        else:
            props = {}

        props['seq_description'] = result['seq_description']

        if props['seq_description'] == 'query':
            props['seq_description'] = 0  # for consistency

        if 'cigar' in result.keys():
            # we did alignment
            props['sequence'] = result['sequence']
            props['score'] = result['score']
            props['cigar'] = result['cigar']

        return {**d, **props}

    lst = [_build_dict(row, result) for result in json_res for row in
           result['results']]

    if lst:
        return pd.DataFrame(lst)
    else:
        # columns may vary on the graph, so not adding column information
        return pd.DataFrame()


def df_from_align_result(json_res):
    # flatten out json result
    lst = [alignment for result in json_res for alignment in
           result['alignments']]

    df = pd.DataFrame(lst,
                      columns=['cigar', 'score', 'sequence', 'seq_description'])

    # for consistency, set seq_description to 0 even if we only queried a single sequence
    df.loc[df['seq_description'] == 'query', 'seq_description'] = 0

    return df
