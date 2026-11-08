# Writing guidance and text checks

Use the [writing principles](https://github.com/ratschlab/writing-rules/blob/fe9defad9a3552bce895c3720243adc61178d6ce/PRINCIPLES.md) for documentation and interface text, and the [engineering communication guide](https://github.com/ratschlab/writing-rules/blob/fe9defad9a3552bce895c3720243adc61178d6ce/ENGINEERING.md) for pull requests, reviews and commit messages.
The `tools/writing-rules` submodule pins these guides and the checker together.

## Run the check

Use Python 3.12. Download the [Vale 3.26.0 archive for your platform](https://github.com/vale-cli/vale/releases/tag/v3.26.0), extract it, and put the `vale` executable in a directory on your `PATH`.
Run `vale --version` to confirm that it reports 3.26.0 before continuing.
CI verifies the Linux x86_64 Vale archive checksum. From the repository root:

```sh
git submodule update --init -- tools/writing-rules
python3.12 -m venv .writing/venv
.writing/venv/bin/pip install ./tools/writing-rules
.writing/venv/bin/writing-rules export-package .writing/WritingRules.zip
vale sync --no-global
PATH="$PWD/.writing/venv/bin:$PATH" .writing/venv/bin/writing-rules check
```

After pulling a new submodule revision, repeat the submodule update and package
installation. Move the old `.writing/WritingRules.zip` aside before exporting
its replacement, then run `vale sync --no-global` again. Export refuses to
overwrite an existing archive. Keep the virtual environment on PATH so Vale can
find Docutils when parsing reStructuredText.

## Interpret the result

The file scope is `writing-rules.json`, and terminology is configured in
`styles/MetaGraph/R8.yml`. The baseline in `writing-rules-baseline.json` records
error counts. The check rejects new terminology errors and malformed exceptions.
Sentence-length, passive-voice and vague-error findings are prompts for review;
judge them in context and preserve the intended meaning.

`writing-rules check --report` reports findings without enforcing the baseline.
`writing-rules check --update-baseline` records reductions in error counts and
refuses increases.

Prefer correcting terminology or marking a literal identifier as code. When a
term must remain in prose, explain the exception in a same-line or preceding
comment, for example `text-rules-allow: R8 kmer the published identifier`.
The rule ID, exact term and reason are required. Native Vale suppression
directives are rejected. This repository currently needs no prose exceptions.

C++ findings point to the message's starting line. Syntax-recovery warnings
identify source that needs manual inspection; the extractor does not perform
data-flow analysis. The current warnings concern config.cpp, query.cpp,
query.hpp, stats.cpp and transform_annotation.cpp. A passing check does not
establish factual correctness, complete extraction coverage, accessibility or
layout quality.
