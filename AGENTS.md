# AGENTS.md

## Writing for people

From a plain clone, run `git submodule update --init -- tools/writing-rules`.
Before changing documentation or user-facing text, read
[`tools/writing-rules/guides/PRINCIPLES.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/guides/PRINCIPLES.md). The accompanying
[`tools/writing-rules/sources/README.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/sources/README.md) explains the basis and limitations of the guidance. The submodule pins the guide and checker together.
For text that explains a concept or algorithm, such as a design note or a
conceptual tutorial, also read
[`tools/writing-rules/guides/TEACHING.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/guides/TEACHING.md).

Write for the reader's task. Preserve facts, uncertainty and necessary technical
detail. Check defaults, flags, file names and guarantees against their owning
code or configuration. Use the actual CLI and API names in examples.

Use “k-mer” consistently in prose. The local vocabulary is
`styles/MetaGraph/R8.yml`; code identifiers and literal command syntax retain
their spelling. Preserve uncertainty and necessary technical detail when
considering a style suggestion.

The selected scope is README, documentation and CLI help/errors, as listed in
`writing-rules.json`. Comments and debug logs are outside the selected C++
message extraction. The same principles guide PR descriptions and reviews,
which are not automatically linted. Describe the concrete change and useful
evidence concisely; report material review findings without a mandatory verdict
phrase or a fixed word count.

Before drafting an issue, PR description, review or commit message, read the relevant
section of [`tools/writing-rules/guides/ENGINEERING.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/guides/ENGINEERING.md).
Its guidance also covers replies to users and release/migration notes.
MetaGraph adopts the Must, Should and Nit scale that guide defines for review
findings. Automated reviewers
such as cubic keep their own P0–P3 priority and add no second label. Read P0 and
P1 as Must, P2 as Should and P3 as Nit.
Before submission, check the receiving project's contribution guidelines,
templates and required disclosures.

Run the [text check](docs/writing/README.md) after editing scoped text.
It blocks increases in terminology errors and malformed exceptions. Sentence
length, passive voice and vague-error phrases are review prompts. Parser
recovery warnings identify coverage that needs inspection. A passing check does
not establish factual correctness or complete extraction coverage.

## Adding features

Before adding a feature, option, flag, endpoint or API field, or reviewing a
change that adds one, read
[`tools/writing-rules/guides/STEWARDSHIP.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/guides/STEWARDSHIP.md).
Name who needs it, what happens without it and the smallest version that serves
the need. Record the decision in the PR.

## Code written with LLM agents

Before submitting or reviewing code written with an LLM agent, run the checklists in
[`tools/writing-rules/guides/LLM_CODE.md`](https://github.com/ratschlab/writing-rules/blob/376b4f302cf4696e0b3a7184efc429da8438be1c/guides/LLM_CODE.md).
For every new function, type or test helper, first search master and the
dependencies for one that already does the job, and call it instead.
