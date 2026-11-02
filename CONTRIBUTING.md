# Contributing

Open a pull request against `master`. The developer notes in the README cover Docker builds and releases.

## Writing user-facing text

These rules apply to text a user reads: `README.md`, `metagraph/docs/source/`, the `--help` text in
`metagraph/src/cli/config/config.cpp`, and the errors and warnings the CLI prints. They don't apply to code
comments, debug logs or commit messages.

1. One idea per sentence. No semicolons. Split the sentence or use a list.
2. Keep sentences short: about 20 words on average, and none over 30.
3. Put the outcome first. The first words say what happened or what the reader must do.
4. An error says what happened, why, and what to do next, in that order. Put codes and IDs after the sentence. Don't
   write "sorry", "oops" or "something went wrong" on its own.
5. Use the reader's words. Domain terms such as *k*-mer, de Bruijn graph or annotation are fine. Define an internal
   name (a format such as `RowDiff<Multi-BRWT>`, a class name) the first time, or leave it to the docs.
6. Use active voice and name the actor: you, MetaGraph, the server, the workflow.
7. Write numbers as digits with units, such as "1 hour" or "5 GB". Don't write "may", "soon" or "a while" when the
   code knows the answer.
8. Use one term per thing and one verb per action. The list is in `scripts/text-rules.json`. Write "k-mer", not
   "kmer". A flag or default in the docs matches `--help` word for word.
9. Don't write "we" unless the text says who "we" is.
10. Every number, limit, default or guarantee in the text matches the code that sets it. A change to the value
    changes the text in the same commit.
11. A `--help` line fits on 1 line and ends with the default in brackets, such as `[0.7]`.

Check: `python3 scripts/check_text_rules.py`. It reads `scripts/text-rules.json` and compares each file with
`scripts/text-rules-baseline.json`. It rejects new hits for R1, R4, R5, R7, R8 and R9. R2 and R6 are advisory
in `--report`. Reviewers check R3, R10 and R11, and the parts of the other rules that need judgment.
Existing baseline hits are allowed. A file that loses a checked hit fails too, so lower its baseline in the same
commit with `--update-baseline`. That command refuses increases. Changes to the scanner's scope need a reviewed
baseline migration. Mark a deliberate exception with a
comment `text-rules-allow: R1 <reason>` on that line or the line before.
A review of a change to user-facing text ends with a line "Text rules: pass" or the list of breaks.
