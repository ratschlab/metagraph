# Writing guidance and checks

`PRINCIPLES.md`, `SOURCES.md`, `ENGINEERING.md`, `LICENSE` and `AUTHORS` are a tracked snapshot from
[ratschlab/writing-rules](https://github.com/ratschlab/writing-rules), at the
commit in `.writing-rules-revision`. Keep the snapshot and tool pin together.
The license here covers the shared guidance; MetaGraph's license is unchanged.

## Local check

Use Python 3.12 and Vale 3.26.0. Install Vale from its official release;
CI verifies the Linux x86_64 archive checksum. From the repository root:

```sh
python3 -m venv .writing/venv
.writing/venv/bin/pip install "git+https://github.com/ratschlab/writing-rules.git@$(cat .writing-rules-revision)"
.writing/venv/bin/writing-rules export-package .writing/WritingRules.zip
vale sync --no-global
PATH="$PWD/.writing/venv/bin:$PATH" .writing/venv/bin/writing-rules check
```

Export creates a new archive and refuses to overwrite one. For an upgrade,
move the previous archive aside before exporting the new one. Docutils from
the environment must be on PATH for RST parsing.

The file scope is `writing-rules.json`, local terminology is in
`styles/MetaGraph/R8.yml`, and error counts are in `writing-rules-baseline.json`.
The reviewed migration has no terminology errors. Earlier semicolon and modal
word counts are not carried forward because those bans are no longer policy.
R2/R6 suggestions and R4 warnings do not block CI.

A deliberate terminology exception uses a same-line or preceding comment:
`text-rules-allow: R8 kmer the published identifier`. The reason is required and
the term must match the finding. Native Vale suppression directives are rejected.
`writing-rules check --report` reports findings without blocking on count increases;
`writing-rules check --update-baseline` records reductions only.

C++ findings point to the message's starting line. Syntax-recovery warnings
remain visible for config.cpp, query.cpp, query.hpp, stats.cpp and
transform_annotation.cpp in the migration corpus. Inspect affected changes;
the extractor does not perform data-flow analysis. Current checks do not test
web UI layout, accessibility or agent output quality.
