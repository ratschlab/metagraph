# -*- coding: utf-8 -*-

import multiprocessing
from typing import Dict, Tuple, List, Iterable, Union, Any
from concurrent.futures import ThreadPoolExecutor, Future, wait

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
               with_signature: bool = False,
               abundance_sum: bool = False,
               query_counts: bool = False,
               query_coords: bool = False,
               align: bool = False,
               **align_params) -> Tuple[JsonDict, str]:
        """See parameters for alignment `align_params` in align()"""

        if discovery_threshold < 0.0 or discovery_threshold > 1.0:
            raise ValueError(
                f"discovery_threshold should be between 0 and 1 inclusive. Got {discovery_threshold}")

        if align:
            # Warn if number of alignments is specified > 1
            if 'max_alternative_alignments' in align_params and \
                    align_params['max_alternative_alignments'] > 1:
                align_params['max_alternative_alignments'] = 1
                warnings.warn(f"Requested max alternative alignments > 1, treating as 1.",
                              RuntimeWarning)
            alignments = self.align(sequence, **align_params)

            # For each input sequence, retrieve the best (first) aligned sequence
            aligned_sequences = []

            for alignment in alignments:
                if len(alignment['alignments']) == 0:
                    # Produce an empty sequence if there are no alignments
                    aligned_sequences.append('')
                else:
                    aligned_sequences.append(alignment['alignments'][0]['sequence'])

            sequence = aligned_sequences

        param_dict = {"count_labels": True,
                      "discovery_fraction": discovery_threshold,
                      "num_labels": top_labels,
                      "with_signature": with_signature,
                      "abundance_sum": abundance_sum,
                      "query_counts": query_counts,
                      "query_coords": query_coords}

        search_results = self._json_seq_query(sequence, param_dict, "search")

        if align:
            assert len(alignments) == len(search_results)

            # Zip best alignment results
            for alignment, search_result in zip(alignments, search_results):
                if 'alignments' in alignment and len(alignment['alignments']) > 0:
                    search_result['best_alignment'] = alignment['alignments'][0]
                else:
                    search_result['best_alignment'] = {}

        else:
            search_results = self._json_seq_query(sequence, param_dict, "search")

        return search_results

    def align(self, sequence: Union[str, Iterable[str]],
              min_exact_match: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Tuple[JsonDict, str]:
        if min_exact_match < 0.0 or min_exact_match > 1.0:
            raise ValueError(
                f"min_exact_match should be between 0 and 1 inclusive. Got {min_exact_match}")

        if max_num_nodes_per_seq_char < 0:
            warnings.warn("max_num_nodes_per_seq_char < 0, treating as infinite", RuntimeWarning)

        params = {'max_alternative_alignments': max_alternative_alignments,
                  'max_num_nodes_per_seq_char': max_num_nodes_per_seq_char,
                  'min_exact_match': min_exact_match}

        return self._json_seq_query(sequence, params, "align")

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
            raise RuntimeError(
                f"Error while calling the server API. {str(ret.status_code)}: {ret.text}")

        if not ret.ok:
            error_msg = json_obj['error'] if 'error' in json_obj.keys() else str(json_obj)
            raise RuntimeError(
                f"Error while calling the server API. {str(ret.status_code)}: {error_msg}")

        return json_obj

    # noinspection PyTypeChecker
    def column_labels(self) -> Tuple[JsonStrList, str]:
        """
        :returns:   All column labels
        :rtype:     JSON list
        """
        return self._do_request("column_labels", {}, post_req=False)

    def stats(self) -> Tuple[dict, str]:
        return self._do_request("stats", {}, post_req=False)

    def ready(self) -> bool:
        try:
            self.stats()
            return True
        except RuntimeError as e:
            if "503: Server is currently initializing" in str(e):
                return False
            raise e


class GraphClient:
    def __init__(self, host: str, port: int, name: str = None, api_path: str = None):
        self._json_client = GraphClientJson(host, port, name, api_path=api_path)
        self.name = self._json_client.name

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               with_signature: bool = False,
               abundance_sum: bool = False,
               query_counts: bool = False,
               query_coords: bool = False,
               align: bool = False,
               **align_params) -> pd.DataFrame:
        """
        :param      sequence:             The query sequence
        :type       sequence:             Union[str, Iterable[str]]
        :param      top_labels:           The maximum number of matched labels to retrieve [default: 100]
        :type       top_labels:           int
        :param      discovery_threshold:  The minimum fraction (between 0.0 and 1.0) of k-mers from the query required to match a label (occur in a sample) in order for that label to show up in the result [default: 0.0]
        :type       discovery_threshold:  float
        :param      with_signature:       Return the signature of k-mer matches
        :type       with_signature:       bool
        :param      abundance_sum:        Compute the sum of abundances for all k-mers matched
        :type       abundance_sum:        bool
        :param      query_counts:         Query k-mer counts
        :type       query_counts:         bool
        :param      query_coords:         Query k-mer coordinates
        :type       query_coords:         bool
        :param      align:                Align the query sequence to the joint graph and query labels for that alignment instead of the original sequence
        :type       align:                bool
        :param      align_params:         The parameters for alignment (see method align())
        :type       align_params:         dictionary

        :return:    A data frame with query results
        :rtype:     pandas.DataFrame
        """

        json_obj = self._json_client.search(sequence, top_labels,
                                            discovery_threshold, with_signature,
                                            abundance_sum, query_counts, query_coords,
                                            align, **align_params)

        return helpers.df_from_search_result(json_obj)

    def align(self, sequence: Union[str, Iterable[str]],
              min_exact_match: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> pd.DataFrame:
        """
        Align sequence(s) to the joint graph

        :param      sequence:                    The query sequence
        :type       sequence:                    Union[str, Iterable[str]]
        :param      min_exact_match:             The minimum fraction (between 0.0 and 1.0) of matching nucleotides required to align sequence [default: 0]
        :type       min_exact_match:             float
        :param      max_alternative_alignments:  The number of different alignments to return [default: 1]
        :type       max_alternative_alignments:  int
        :param      max_num_nodes_per_seq_char:  The maximum number of nodes to consider per sequence character during extension [default: 10.0]
        :type       max_num_nodes_per_seq_char:  float

        :returns:   A data frame with alignments
        :rtype:     pandas.DataFrame
        """
        json_obj = self._json_client.align(sequence, min_exact_match,
                                           max_alternative_alignments,
                                           max_num_nodes_per_seq_char)

        return helpers.df_from_align_result(json_obj)

    def column_labels(self) -> List[str]:
        return self._json_client.column_labels()

    def ready(self) -> bool:
        return self._json_client.ready()


class MultiGraphClient:
    """
    A version of the graph client that can add multiple graph servers and also supports
    asynchronous querying using a thread pool.
    """
    def __init__(self):
        """ Create an instance of MultiGraphClient. """
        self.graphs = {}

    def add_graph(self, host: str, port: int, name: str = None, api_path: str = None) -> None:
        """ Adds graph client to list of graphs to query on request """

        graph_client = GraphClient(host, port, name, api_path=api_path)
        self.graphs[graph_client.name] = graph_client

    def list_graphs(self) -> Dict[str, Tuple[str, int]]:
        """ List the details of the currently added graph clients """

        return {lbl: (inst._json_client.host, inst._json_client.port) for (lbl, inst) in
                self.graphs.items()}

    def search(self, sequence: Union[str, Iterable[str]],
               parallel=True,
               top_labels: int = DEFAULT_TOP_LABELS,
               discovery_threshold: float = DEFAULT_DISCOVERY_THRESHOLD,
               with_signature: bool = False,
               abundance_sum: bool = False,
               query_counts: bool = False,
               query_coords: bool = False,
               align: bool = False,
               **align_params) -> Dict[str, Union[pd.DataFrame, Future]]:
        """
        Make search requests with each graph client.

        :param parallel: perform requests in parallel, and instead return Future instances
        :param align_params: see parameters for alignment in align()
        :return: Dict mapping graph names to pd.DataFrame results
        """

        if not parallel:
            result = {}

            # Do this iteratively and return once done
            for name, graph_client in self.graphs.items():
                result[name] = graph_client.search(sequence, top_labels,
                                                   discovery_threshold, with_signature,
                                                   abundance_sum, query_counts, query_coords,
                                                   align, **align_params)

            return result

        # Otherwise, let's query in parallel with a multiprocessing pool
        num_processes = min(multiprocessing.cpu_count(), len(self.graphs))

        futures = {}
        executor = ThreadPoolExecutor(max_workers=num_processes)

        # Populate async results dict with concurrent.futures.Future instances
        for name, graph_client in self.graphs.items():
            futures[name] = executor.submit(graph_client.search, sequence,
                                            top_labels, discovery_threshold, with_signature,
                                            abundance_sum, query_counts, query_coords,
                                            align, **align_params)

        print(f'Made {len(self.graphs)} requests with {num_processes} threads...')

        # Shutdown executor but do not stop futures
        executor.shutdown(wait=False)
        return futures

    def align(self, sequence: Union[str, Iterable[str]],
              parallel=True,
              min_exact_match: float = DEFAULT_DISCOVERY_THRESHOLD,
              max_alternative_alignments: int = 1,
              max_num_nodes_per_seq_char: float = DEFAULT_NUM_NODES_PER_SEQ_CHAR) -> Dict[
        str, Union[pd.DataFrame, Future]]:
        """
        Make align requests with each graph client.

        :param parallel: perform requests in parallel, and instead return Future instances
        :return: Dict mapping graph names to pd.DataFrame results
        """

        if not parallel:
            result = {}

            # Do this iteratively and return once done
            for name, graph_client in self.graphs.items():
                result[name] = graph_client.align(sequence, min_exact_match,
                                                  max_alternative_alignments,
                                                  max_num_nodes_per_seq_char)

            return result

        # Otherwise, let's query in parallel with a multiprocessing pool
        num_processes = min(multiprocessing.cpu_count(), len(self.graphs))

        futures = {}
        executor = ThreadPoolExecutor(max_workers=num_processes)

        # Populate async results dict with concurrent.futures.Future instances
        for name, graph_client in self.graphs.items():
            futures[name] = executor.submit(graph_client.align, sequence, min_exact_match,
                                            max_alternative_alignments,
                                            max_num_nodes_per_seq_char)

        print(f'Made {len(self.graphs)} requests with {num_processes} threads...')

        # Shutdown executor but do not stop futures
        executor.shutdown(wait=False)
        return futures

    def column_labels(self) -> Dict[str, List[str]]:
        ret = {}
        for name, graph_client in self.graphs.items():
            ret[name] = graph_client.column_labels()

        return ret

    def wait_for_result(futures: Dict[str, Future]) -> Dict[str, Union[pd.DataFrame, Exception]]:
        """
        Wait for all futures to finish running and return a proper result dict.
        If any of the calls resulted in an exception, return that in the result instead.

        :param futures: Dictionary mapping graph names to Future instances
        :return Dictionary mapping graph names to result DataFrames (or Execeptions)
        """
        wait(futures.values())

        result = {}
        for name, fs in futures.items():
            result[name] = fs.exception() if fs.exception() else fs.result()

        return result
