# -*- coding: utf-8 -*-

import json
from typing import Dict, Tuple, List

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

    # TODO: discovery threshoold unit?
    def search(self, sequence: str, discovery_threshold: float = 1.0) -> Dict[
        str, Tuple[pd.DataFrame, str]]:

        json_res = self.search_json(sequence, discovery_threshold)

        def build_dict(row):
            d = dict(row)
            props = d.pop('properties')
            return {**d, **props}

        def build_df_from_json(j):
            return pd.DataFrame([build_dict(r) for r in j['results']])

        return {l: (build_df_from_json(j) if j else pd.DataFrame()) for (l, (j, e))
                in json_res.items()}

    def search_json(self, sequence: str, discovery_threshold: float = 1.0) -> \
    Dict[str, Tuple[str, str]]:
        if not self.graphs:
            raise ValueError("No graphs registered")  # TODO: better error

        results = {}

        for label, (host, port) in self.graphs.items():
            payload = json.dumps({
                "FASTA": "\n".join([">query",
                                    sequence,
                                    ]),
                "count_labels": True,
                "discovery_fraction": discovery_threshold / 100,
                "num_labels": 10000
            })

            ret = requests.post(url=f'http://{host}:{port}/search', data=payload)

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
