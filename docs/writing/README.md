# Writing guidance and checks

The [writing-rules submodule](../../tools/writing-rules) contains the guide and
checker at the revision selected by MetaGraph. Initialize it once after cloning,
and repeat the update command after pulling a new submodule revision.

## Local check

Use Python 3.12 and Vale 3.26.0. Install Vale from its official release;
CI verifies the Linux x86_64 archive checksum. From the repository root:

```sh
git submodule update --init -- tools/writing-rules
python3 -m venv .writing/venv
.writing/venv/bin/pip install ./tools/writing-rules
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
