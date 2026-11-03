# AGENTS.md

Rules for coding agents working in this repository.

## Writing for people

These rules apply to everything a person reads in this repository:

- `README.md` and the docs in `metagraph/docs/source/`.
- The Python API README at `metagraph/api/python/README.rst`.
- The `--help` text in `metagraph/src/cli/config/config.cpp`.
- The errors and warnings the CLI prints.
- PR descriptions, review comments and commit messages.

They don't apply to code comments or debug logs. The numbers follow the shared list R1 to R28. A missing number is a
rule that doesn't apply here.

- R1. One idea per sentence. No semicolons, no em dash between clauses, and at most 1 colon per sentence. Split the
  sentence or use a list.
- R2. Short sentences: 20 words or fewer on average, none over 30.
- R3. Outcome first. The first 2 words carry the outcome or the action.
- R4. An error says what happened, why, and what to do next, in that order. The next step names who or what: a flag,
  a file, a page, a contact. Codes and IDs go after the sentence. No "sorry", "oops" or a bare "something went
  wrong".
- R5. The reader's words. Domain terms such as *k*-mer, de Bruijn graph or annotation are fine. Define an internal
  name (a format such as `RowDiff<Multi-BRWT>`, a class name) the first time, or leave it to the docs.
- R6. Active voice with a named actor: you, MetaGraph, the server, the workflow. Passive voice is fine when no actor
  matters, as in "Nothing was changed."
- R7. Digits and units for every number and time: "1 hour", "5 GB". Write "may" only for permission. No "soon" or "a
  while" when the code knows the answer.
- R8. One term per thing and one verb per action. The term list in `scripts/text-rules.json` is closed: a new term
  goes into it in the same commit. Write "k-mer", not "kmer". A flag or default in the docs matches `--help` word
  for word.
- R9. No "we" unless the text says who "we" is.
- R10. Every number, limit, default or guarantee matches the code or config that sets it. The commit that adds a
  value to the text names where the code sets it. A change to the value changes the text in the same commit.
- R12. A `--help` option fits on 1 line. An option with a default ends with it in brackets, such as `[0.7]`.
- R17. Paragraphs: at most 5 sentences in docs and 3 sentences in an error. Three or more parallel items become a
  list.
- R18. Headings are statements that name the topic, with the key word first.
- R19. Details last. IDs, codes, raw output and internals go after what most readers need, in a Details part.
- R23. A link names its destination, never "here" or "click here".

What a PR, a review and a commit look like:

- **Commit message.** The subject says what changed, in plain words, in 72 characters or fewer. The body says why,
  then the numbers (before and after, with the command), then the risk.
- **PR description.** At most about 15 lines. The first line is the outcome for a user. Then 2 to 4 bullets: what
  changed and why. Then the numbers with the command that produced them, the risk, and how it was checked (`Tests:
  CI`, since tests run on CI only). Say what is not done. Don't list open review points.
- **Review comment.** One finding per comment. Start with the severity (Must, Should or Nit) and the outcome: what
  breaks, for whom. Then the evidence (`file:line`, input, output) and a concrete fix. End the review with a line
  `Text rules: pass` or the rule numbers with `file:line`.

`python3 scripts/check_text_rules.py` compares the files in `scripts/text-rules.json` with
`scripts/text-rules-baseline.json`. It fails only when a file gains one of these hits:

- R1: a semicolon in prose.
- R4: "sorry", "oops", "whoops", "something went wrong", "an error occurred" or "unknown error".
- R5: a word from the jargon list, in prose.
- R7: "may", "might", "soon", "shortly", "a while" or "eventually".
  It also detects selected number words before selected units, ignoring case.
  Numbers: two, three, four, five, six, seven, eight, nine, ten, twelve, twenty, thirty and sixty.
  Units: second, minute, hour, day, week, month, thread, file, step and result, including their plurals.
  It also recognizes times, GB, MB and TB. Other number words and units need manual review.
- R8: a banned term from the term list.
- R9: "we", "our" or "us" where the scope bans them.

It reports R2 (sentences over 30 words) and R6 (passive voice) with `--report` and never fails on them. Every other
rule, and em dashes and colons for R1, is for the reviewer. Text outside these files, such as commit messages and PR
descriptions, is not linted.

A file with fewer hits than its baseline passes with a note. Lower its baseline with `--update-baseline`, which
refuses increases. A deliberate exception takes a comment on that line or the line before: `text-rules-allow: <rule>
<term> <reason>`, such as `text-rules-allow: R8 kmer the flag name`. A reason is required.
The term must occur in the text with the same word boundaries used by the term check.

Reviewers use this list.
