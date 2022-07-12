import pandas as pd


def df_from_search_result(json_res):
    def _build_dict(row, result):
        d = dict(row)
        # `properties` (optional): dictionary with metadata about the sample
        if 'properties' in d.keys():
            props = d.pop('properties')
        else:
            props = {}

        props['seq_description'] = result['seq_description']

        return {**d, **props}

    lst = [_build_dict(row, result) for result in json_res for row in
           result['results']]

    # columns of the table may vary on the graph, and inferred automatically
    return pd.DataFrame(lst)


def df_from_align_result(json_res):
    # flatten out json result
    lst = [(alignment['cigar'],
            alignment['score'],
            alignment['max_score'],
            alignment['sequence'],
            alignment['orientation'],
            result['seq_description'])
           for result in json_res for alignment in result['alignments']]

    df = pd.DataFrame(lst,
                      columns=['cigar', 'score', 'max_score', 'sequence', 'orientation', 'seq_description'])

    return df
