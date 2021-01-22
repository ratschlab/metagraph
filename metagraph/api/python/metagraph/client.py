# -*- coding: utf-8 -*-

from typing import Dict, Tuple, List, Iterable, Union, Any

import pandas as pd
import requests
import warnings

from metagraph import helpers

"""Metagraph client."""

DEFAULT_TOP_LABELS = 100
DEFAULT_DISCOVERY_THRESHOLD = 0
DEFAULT_NUM_NODES_PER_SEQ_CHAR = 10.0

JsonDict = Dict[str, Any]
JsonStrList = List[str]


class GraphClientJson:
    """
    Relatively low level version of the client API. Client returning results
    from the server as json objects. If there was an error,
    returning error message in the second element of the tuple returned.
    """

    def __init__(self, host: str, port: int, name: str = None, api_path: str = None):
        self.host = host
        self.port = port

        self.server = f"http://{self.host}:{self.port}"
        if api_path:
            self.server = f"{self.server}/{api_path.lstrip('/')}"

        self.name = name if name else self.server

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               **align_params) -> Tuple[JsonDict, str]:
        """See parameters for alignment `align_params` in align()"""

        if discovery_threshold < 0.0 or discovery_threshold > 1.0:
            raise ValueError(
                f"discovery_threshold should be between 0 and 1 inclusive. Got {discovery_threshold}")

        if align:
            json_obj, err = self.align(sequence, discovery_threshold, **align_params)

            if err:
                raise RuntimeError(f"Error while calling the server API {str(err)}")

            def to_fasta(df):
                fasta = []
                for i in range(df.shape[0]):
                    fasta.append(f">{df.loc[i, 'seq_description']}\n{df.loc[i, 'sequence']}")
                return '\n'.join(fasta)

            sequence = to_fasta(helpers.df_from_align_result(json_obj))

        param_dict = {"count_labels": True,
                      "discovery_fraction": discovery_threshold,
                      "num_labels": top_labels}

        return self._json_seq_query(sequence, param_dict, "search")

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Tuple[JsonDict, str]:
        if discovery_threshold < 0.0 or discovery_threshold > 1.0:
            raise ValueError(
                f"discovery_threshold should be between 0 and 1 inclusive. Got {discovery_threshold}")

        if max_num_nodes_per_seq_char < 0:
            warnings.warn("max_num_nodes_per_seq_char < 0, treating as infinite", RuntimeWarning)

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
            sequence = [sequence]

        fasta_str = '\n'.join(
            [f">{i}\n{seq}" for (i, seq) in enumerate(sequence)])

        payload_dict = {"FASTA": fasta_str}
        payload_dict.update(param_dict)
        payload = payload_dict

        return self._do_request(endpoint, payload)

    def _do_request(self, endpoint, payload, post_req=True) -> Tuple[JsonDict, str]:
        url = f'{self.server}/{endpoint}'
        if post_req:
            ret = requests.post(url=url, json=payload)
        else:
            ret = requests.get(url=url)

        try:
            json_obj = ret.json()
        except:
            return {}, str(ret.status_code) + " " + ret.text

        if not ret.ok:
            error_msg = json_obj['error'] if 'error' in json_obj.keys() else str(json_obj)
            return {}, str(ret.status_code) + " " + error_msg

        return json_obj, ""

    def stats(self) -> Tuple[dict, str]:
        return self._do_request("stats", {}, post_req=False)

    def ready(self) -> bool:
        result = self.stats()
        if result[0]:
            return True
        if result[1].startswith("503 Server is currently initializing"):
            return False
        raise RuntimeError(f"Error. Server response: {result[1]}")


class GraphClient:
    def __init__(self, host: str, port: int, name: str = None, api_path: str = None):
        self._json_client = GraphClientJson(host, port, name, api_path=api_path)
        self.name = self._json_client.name

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               **align_params) -> pd.DataFrame:
        """See parameters for alignment `align_params` in align()"""

        json_obj, err = self._json_client.search(sequence, top_labels,
                                                 discovery_threshold,
                                                 align, **align_params)

        if err:
            raise RuntimeError(
                f"Error while calling the server API {str(err)}")

        return helpers.df_from_search_result(json_obj)

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> pd.DataFrame:
        json_obj, err = self._json_client.align(sequence, discovery_threshold,
                                                max_alternative_alignments,
                                                max_num_nodes_per_seq_char)

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")

        return helpers.df_from_align_result(json_obj)

    def column_labels(self) -> List[str]:
        json_obj, err = self._json_client.column_labels()

        if err:
            raise RuntimeError(f"Error while calling the server API {str(err)}")
        return json_obj

    def ready(self) -> bool:
        return self._json_client.ready()


class MultiGraphClient:
    # TODO: make things asynchronously. this should be the added value of this class
    def __init__(self):
        self.graphs = {}

    def add_graph(self, host: str, port: int, name: str = None, api_path: str = None) -> None:
        graph_client = GraphClient(host, port, name, api_path=api_path)
        self.graphs[graph_client.name] = graph_client

    def list_graphs(self) -> Dict[str, Tuple[str, int]]:
        return {lbl: (inst.host, inst.port) for (lbl, inst) in
                self.graphs.items()}

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               align: bool = False,
               **align_params) -> Dict[str, pd.DataFrame]:
        """See parameters for alignment `align_params` in align()"""

        result = {}
        for name, graph_client in self.graphs.items():
            result[name] = graph_client.search(sequence, top_labels,
                                                discovery_threshold,
                                                align, **align_params)

        return result

    def align(self, sequence: Union[str, Iterable[str]],
              discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Dict[
        str, pd.DataFrame]:
        result = {}
        for name, graph_client in self.graphs.items():
            # TODO: do this async
            result[name] = graph_client.align(sequence, discovery_threshold,
                                               max_alternative_alignments,
                                               max_num_nodes_per_seq_char)

        return result

    def column_labels(self) -> Dict[str, List[str]]:
        ret = {}
        for name, graph_client in self.graphs.items():
            ret[name] = graph_client.column_labels()

        return ret
