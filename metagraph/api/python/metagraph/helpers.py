import pandas as pd


def df_from_search_result(json_res):
    def _build_dict(row, result):
        d = dict(row)
        if 'properties' in d.keys():
            props = d.pop('properties')
        else:
            props = {}

        props['seq_description'] = result['seq_description']

        # TODO: remove
        if 'cigar' in result.keys():
            # we did alignment
            props['sequence'] = result['sequence']
            props['score'] = result['score']
            props['cigar'] = result['cigar']

        return {**d, **props}

    lst = [_build_dict(row, result) for result in json_res for row in
           result['results']]

    # columns of the table may vary on the graph, and inferred automatically
    return pd.DataFrame(lst)


def df_from_align_result(json_res):
    # flatten out json result
    lst = [(alignment['cigar'],
            alignment['score'],
            alignment['sequence'],
            result['seq_description'])
           for result in json_res for alignment in result['alignments']]

    df = pd.DataFrame(lst,
                      columns=['cigar', 'score', 'sequence', 'seq_description'])

    return df
