# -*- coding: utf-8 -*-

import json
from typing import Dict, Tuple, List, Iterable, Union

import pandas as pd
import requests

"""Metagraph client."""


class Client:
    def __init__(self):
        self.graphs = {}

    def add_graph(self, host: str, port: int, label: str = None) -> None:
        if not label:
            label = f"{host}:{port}"
        self.graphs[label] = (host, port)

    def list_graphs(self) -> Dict[str, Tuple[str, str]]:
        return self.graphs

    def search(self, sequence: Union[str, Iterable[str]],
               top_labels: int = 10000, discovery_threshold: float = 1.0, align: bool = False) -> \
            Dict[str, Tuple[pd.DataFrame, str]]:

        json_res = self.search_json(sequence, top_labels, discovery_threshold, align)

        # def build_dict(row):
        #     d = dict(row)
        #     props = d.pop('properties')
        #     return {**d, **props}

        #def build_df_from_json(j):
        #    return pd.DataFrame([build_dict(r) for r in j['results']])

        def _build_df_per_result(res):
            df = pd.DataFrame(res['results'])

            if not isinstance(sequence, str):
                # only add sequence description if several queries are being made
                df['seq_description'] = res['seq_description']

            if align:
                df['sequence'] = res['sequence']
                df['score'] = res['score']

            return df

        def build_df_from_json(j):
            return pd.concat(_build_df_per_result(query_res) for query_res in j)

        return {l: (build_df_from_json(j) if j else pd.DataFrame()) for (l, (j, e))
                in json_res.items()}

    def align(self, sequence: Union[str, Iterable[str]]) -> Dict[
        str, Tuple[pd.DataFrame, str]]:

        json_res = self.align_json(sequence)

        return {l: (pd.DataFrame(j) if j else pd.DataFrame()) for (l, (j, e))
                in json_res.items()}

    def search_json(self, sequence: Union[str, Iterable[str]], top_labels=10000,
                    discovery_threshold: float = 1.0, align: bool = False) -> \
            Dict[str, Tuple[str, str]]:
        param_dict = {"count_labels": True,
                      "discovery_fraction": discovery_threshold / 100,
                      "num_labels": top_labels,
                      "align": align}

        return self._json_query(sequence, param_dict, "search")

    def align_json(self, sequence: Union[str, Iterable[str]]) -> Dict[
        str, Tuple[str, str]]:
        return self._json_query(sequence, {}, "align")

    def _json_query(self, sequence: Union[str, Iterable[str]], param_dict, end_point: str) -> Dict[
        str, Tuple[str, str]]:

        if not self.graphs:
            raise ValueError("No graphs registered")  # TODO: better error

        results = {}

        for label, (host, port) in self.graphs.items():
            if isinstance(sequence, str):
                fasta_str = f">query\n{sequence}"
            else:
                seqs = list(sequence)
                fasta_str = '\n'.join([ f">{i}\n{seqs[i]}" for i in range(0, len(seqs))])

            payload_dict = {"FASTA": fasta_str}
            payload_dict.update(param_dict)
            payload = json.dumps(payload_dict)

            ret = requests.post(url=f'http://{host}:{port}/{end_point}', data=payload)

            if ret.ok:
                results[label] = (ret.json(), None)
            else:
                results[label] = ("", str(ret.status_code) + " " + ret.json()['error'])
        return results

    def column_labels(self) -> Dict[str, List[str]]:
        ret = {}
        for label, (host, port) in self.graphs.items():
            r = requests.get(url=f'http://{host}:{port}/column_labels')

            # TODO: how to handle errors?
            ret[label] = r.json()

        return ret