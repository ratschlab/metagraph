# -*- coding: utf-8 -*-

import json
from typing import Dict, Tuple, List, Iterable, Union, Any

import pandas as pd
import requests

"""Metagraph client."""

DEFAULT_TOP_LABELS = 10000
DEFAULT_DISCOVERY_THRESHOLD = 1.0

JsonDict = Dict[str, Any]
JsonStrList = List[str]


class GraphClientJson:
    """
    Relatively low level version of the client API. Client returning results
    from the server as json objects. If there was an error,
    returning error message in the second element of the tuple returned.
    """

    def __init__(self, host: str, port: int, label: str = None):
        self.host = host
        self.port = port
        self.label = label

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False) -> Tuple[JsonDict, str]:
        if discovery_threshold < 0.0 or discovery_threshold > 1.0:
            raise ValueError(
                f"discovery_threshold should be between 0 and 1 inclusive. Got {discovery_threshold}")

        param_dict = {"count_labels": True,
                      "discovery_fraction": discovery_threshold,
                      "num_labels": top_labels,
                      "align": align}

        return self._json_seq_query(sequence, param_dict, "search")

    def align(self, sequence: Union[str, Iterable[str]]) -> Tuple[JsonDict, str]:
        return self._json_seq_query(sequence, {}, "align")

    # noinspection PyTypeChecker
    def column_labels(self) -> Tuple[JsonStrList, str]:
        return self._do_request("column_labels", {}, False)

    def _json_seq_query(self, sequence: Union[str, Iterable[str]], param_dict,
                        endpoint: str) -> Tuple[JsonDict, str]:
        if isinstance(sequence, str):
            fasta_str = f">query\n{sequence}"
        else:
            seqs = list(sequence)
            fasta_str = '\n'.join(
                [f">{i}\n{seq}" for (i, seq) in enumerate(seqs)])

        payload_dict = {"FASTA": fasta_str}
        payload_dict.update(param_dict)
        payload = json.dumps(payload_dict)

        return self._do_request(endpoint, payload)

    def _do_request(self, endpoint, payload, post_req=True) -> Tuple[JsonDict, str]:
        url = f'http://{self.host}:{self.port}/{endpoint}'
        if post_req:
            ret = requests.post(url=url, data=payload)
        else:
            ret = requests.get(url=url)

        json_obj = ret.json()
        if not ret.ok:
            error_msg = json_obj[
                'error'] if 'error' in json_obj.keys() else str(json_obj)
            return {}, str(ret.status_code) + " " + error_msg

        return json_obj, ""


class GraphClient:
    def __init__(self, host: str, port: int, label: str = None):
        self._json_client = GraphClientJson(host, port)
        self.label = label

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False) -> pd.DataFrame:
        (json_obj, err) = self._json_client.search(sequence, top_labels,
                                                   discovery_threshold, align)

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

            return df

        def build_df_from_json(j):
            return pd.concat(_build_df_per_result(query_res) for query_res in j)

        return build_df_from_json(json_obj)

    def align(self, sequence: Union[str, Iterable[str]]) -> pd.DataFrame:
        json_obj, err = self._json_client.align(sequence)

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")
        return pd.DataFrame(json_obj)

    def column_labels(self) -> List[str]:
        json_obj, err = self._json_client.column_labels()

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")
        return json_obj


class MultiGraphClient:
    # TODO: make things asynchronously. this should be the added value of this class
    def __init__(self):
        self.graphs = {}

    def add_graph(self, host: str, port: int, label: str = None) -> None:
        if not label:
            label = f"{host}:{port}"

        self.graphs[label] = GraphClient(host, port, label)

    def list_graphs(self) -> Dict[str, Tuple[str, int]]:
        return {lbl: (inst.host, inst.port) for (lbl, inst) in
                self.graphs.items()}

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False) -> \
            Dict[str, pd.DataFrame]:

        result = {}
        for label, graph_client in self.graphs.items():
            result[label] = graph_client.search(sequence, top_labels,
                                                discovery_threshold, align)

        return result

    def align(self, sequence: Union[str, Iterable[str]]) -> Dict[
        str, pd.DataFrame]:
        result = {}
        for label, graph_client in self.graphs.items():
            # TODO: do this async
            result[label] = graph_client.align(sequence)

        return result

    def column_labels(self) -> Dict[str, List[str]]:
        ret = {}
        for label, graph_client in self.graphs.items():
            ret[label] = graph_client.column_labels()

        return ret
