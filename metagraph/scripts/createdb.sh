#!/bin/bash

# TODO refactor DB_NAME

createdb succinct_debruijn_graph

psql -d succinct_debruijn_graph -a -f createtables.sql

# dropdb succinct_debruijn_graph
