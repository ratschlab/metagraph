# AGENTS.md

## Writing for people

Before changing documentation or user-facing text, read
[the shared principles](docs/writing/PRINCIPLES.md). The accompanying
[sources](docs/writing/SOURCES.md) distinguish research, standards and style
preferences. `.writing-rules-revision` pins this local snapshot and the checker.

Write for the reader's task. Preserve facts, uncertainty and necessary technical
detail. Check defaults, flags, file names and guarantees against their owning
code or configuration. Use the actual CLI and API names in examples.

Use “k-mer” consistently in prose. The local vocabulary is
`styles/MetaGraph/R8.yml`; code identifiers and literal command syntax retain
their spelling. No blanket ban on semicolons, “may”, “might” or mathematical “we”
applies. Do not change meaning to satisfy a style suggestion.

The selected scope is README, documentation and CLI help/errors, as listed in
`writing-rules.json`. Comments and debug logs are outside the selected C++
message extraction. The same principles guide PR descriptions and reviews,
which are not automatically linted. Describe the concrete change and useful
evidence concisely; report material review findings without a mandatory verdict
phrase or a fixed word count.

Before drafting a PR description, review or commit message, read the relevant
section of [engineering communication](docs/writing/ENGINEERING.md).
Its guidance also covers bug reports and release/migration notes.

Run the [text check](docs/writing/README.md) after editing scoped text.
It blocks increases in terminology errors and malformed exceptions. Sentence
length, passive voice and vague-error phrases are review prompts. Parser
recovery warnings identify coverage that needs inspection. A passing check does
not establish factual correctness or complete extraction coverage.
