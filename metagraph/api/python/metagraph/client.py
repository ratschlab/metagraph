# -*- coding: utf-8 -*-

from typing import Dict, Tuple, List, Iterable, Union, Any

import pandas as pd
import requests

"""Metagraph client."""

DEFAULT_TOP_LABELS = 10000
DEFAULT_DISCOVERY_THRESHOLD = 0.7
DEFAULT_NUM_NODES_PER_SEQ_CHAR = 10.0

JsonDict = Dict[str, Any]
JsonStrList = List[str]


class GraphClientJson:
    """
    Relatively low level version of the client API. Client returning results
    from the server as json objects. If there was an error,
    returning error message in the second element of the tuple returned.
    """

    def __init__(self, host: str, port: int, label: str = None, api_path: str = None):
        self.host = host
        self.port = port
        self.label = label
        self.api_path = api_path

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Tuple[JsonDict, str]:
        if discovery_threshold < 0.0 or discovery_threshold > 1.0:
            raise ValueError(
                f"discovery_threshold should be between 0 and 1 inclusive. Got {discovery_threshold}")

        param_dict = {"count_labels": True,
                      "discovery_fraction": discovery_threshold,
                      "num_labels": top_labels,
                      "align": align,
                      "max_num_nodes_per_seq_char": max_num_nodes_per_seq_char}

        return self._json_seq_query(sequence, param_dict, "search")

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Tuple[JsonDict, str]:
        params = {'max_alternative_alignments': max_alternative_alignments,
                  'max_num_nodes_per_seq_char': max_num_nodes_per_seq_char,
                  'discovery_fraction': discovery_threshold}
        return self._json_seq_query(sequence, params, "align")

    # noinspection PyTypeChecker
    def column_labels(self) -> Tuple[JsonStrList, str]:
        return self._do_request("column_labels", {}, False)

    def _json_seq_query(self, sequence: Union[str, Iterable[str]], param_dict,
                        endpoint: str) -> Tuple[JsonDict, str]:
        if isinstance(sequence, str):
            fasta_str = f">query\n{sequence}"
        else:
            fasta_str = '\n'.join(
                [f">{i}\n{seq}" for (i, seq) in enumerate(sequence)])

        payload_dict = {"FASTA": fasta_str}
        payload_dict.update(param_dict)
        payload = payload_dict

        return self._do_request(endpoint, payload)

    def _do_request(self, endpoint, payload, post_req=True) -> Tuple[JsonDict, str]:
        endpoint_path = endpoint

        if self.api_path:
            endpoint_path = f"{self.api_path.lstrip('/')}/{endpoint}"

        url = f'http://{self.host}:{self.port}/{endpoint_path}'
        if post_req:
            ret = requests.post(url=url, json=payload)
        else:
            ret = requests.get(url=url)

        try:
            json_obj = ret.json()
        except:
            return {}, str(ret.status_code) + " " + str(ret)

        if not ret.ok:
            error_msg = json_obj[
                'error'] if 'error' in json_obj.keys() else str(json_obj)
            return {}, str(ret.status_code) + " " + error_msg

        return json_obj, ""

    def stats(self) -> Tuple[dict, str]:
        return self._do_request("stats", {}, post_req=False)


class GraphClient:
    def __init__(self, host: str, port: int, label: str = None, api_path: str = None):
        self._json_client = GraphClientJson(host, port, api_path=api_path)
        self.label = label

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> pd.DataFrame:
        (json_obj, err) = self._json_client.search(sequence, top_labels,
                                                   discovery_threshold, align,
                                                   max_num_nodes_per_seq_char)

        if err:
            raise RuntimeError(
                f"Error while calling the server API {str(err)}")

        def _build_dict(row):
            d = dict(row)
            if 'properties' in d.keys():
                props = d.pop('properties')
            else:
                props = {}
            return {**d, **props}

        def _build_df_from_json(j):
            return pd.DataFrame([_build_dict(r) for r in j['results']])

        def _build_df_per_result(res):
            df = _build_df_from_json(res)

            if not isinstance(sequence, str):
                # only add sequence description if several queries are being made
                df['seq_description'] = res['seq_description']

            if align:
                df['sequence'] = res['sequence']
                df['score'] = res['score']
                df['cigar'] = res['cigar']

            return df

        def build_df_from_json(j):
            if j:
                return pd.concat(_build_df_per_result(query_res) for query_res in j)
            return pd.DataFrame()

        return build_df_from_json(json_obj)

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> pd.DataFrame:
        json_obj, err = self._json_client.align(sequence, discovery_threshold,
                                                max_alternative_alignments, max_num_nodes_per_seq_char)

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")

        def _df_per_seq_res(seq_res):
            df = pd.DataFrame(seq_res['alignments'])
            if 'cigar' not in df:
                # no results returned. add empty columns
                df['cigar'] = ""
                df['score'] = None
                df['sequence'] = ""
            df['seq_description'] = seq_res['seq_description']
            return df

        if not json_obj:
            return pd.DataFrame({})

        return (pd.concat([ _df_per_seq_res(a) for a in json_obj]).
                 reset_index(drop=True))


    def column_labels(self) -> List[str]:
        json_obj, err = self._json_client.column_labels()

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")
        return json_obj


class MultiGraphClient:
    # TODO: make things asynchronously. this should be the added value of this class
    def __init__(self):
        self.graphs = {}

    def add_graph(self, host: str, port: int, label: str = None, api_path: str = None) -> None:
        if not label:
            label = f"{host}:{port}"

        self.graphs[label] = GraphClient(host, port, label, api_path=api_path)

    def list_graphs(self) -> Dict[str, Tuple[str, int]]:
        return {lbl: (inst.host, inst.port) for (lbl, inst) in
                self.graphs.items()}

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> \
            Dict[str, pd.DataFrame]:

        result = {}
        for label, graph_client in self.graphs.items():
            result[label] = graph_client.search(sequence, top_labels,
                                                discovery_threshold,
                                                align,
                                                max_num_nodes_per_seq_char)

        return result

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Dict[
        str, pd.DataFrame]:
        result = {}
        for label, graph_client in self.graphs.items():
            # TODO: do this async
            result[label] = graph_client.align(sequence, discovery_threshold,
                                               max_alternative_alignments,
                                               max_num_nodes_per_seq_char)

        return result

    def column_labels(self) -> Dict[str, List[str]]:
        ret = {}
        for label, graph_client in self.graphs.items():
            ret[label] = graph_client.column_labels()

        return ret
