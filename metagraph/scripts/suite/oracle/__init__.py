"""Lane O: a small synthetic-graph suite with an independent brute-force oracle.

The oracle (oracle.py) enumerates the complete answer of a /pattern or /resolve request over the
text of a few tiny records (fixtures.py); compare.py holds an engine answer against it, check by
check; mutations.py corrupts answers one property at a time; run_oracle.py builds the indexes with
the engine's CLI, serves them on localhost, runs the generated requests (requests_gen.py) and
writes the report (report.py). See the section "Small synthetic-graph oracle (lane O)" of
the suite README.

Nothing here imports the engine, checks.py, seqtools.py, truth.py or the integration tests: the
oracle's definitions are written from the SPEC's text.
"""
