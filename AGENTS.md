# AGENTS.md

Rules for coding agents working in this repository.

## Writing for people

These rules apply to everything a person reads in this repository: `README.md`, `metagraph/docs/source/`, the
`--help` text in `metagraph/src/cli/config/config.cpp`, the errors and warnings the CLI prints, and the PR
descriptions, review comments and commit messages of this repository. They don't apply to code comments or debug
logs. The numbers follow the shared list R1 to R28, so a missing number is a rule that doesn't apply here.

- R1. One idea per sentence. No semicolons, no em dash between clauses, and at most 1 colon per sentence. Split the
  sentence or use a list.
- R2. Short sentences: 20 words or fewer on average, none over 30.
- R3. Outcome first. The first 2 words carry the outcome or the action.
- R4. An error says what happened, why, and what to do next, in that order. The next step names who or what: a flag,
  a file, a page, a contact. Codes and IDs go after the sentence. No "sorry", "oops" or a bare "something went
  wrong".
- R5. The reader's words. Domain terms such as *k*-mer, de Bruijn graph or annotation are fine. Define an internal
  name (a format such as `RowDiff<Multi-BRWT>`, a class name) the first time, or leave it to the docs.
- R6. Active voice with a named actor: you, MetaGraph, the server, the workflow. A short outcome with no actor, such
  as "Nothing was changed.", can stay passive.
- R7. Digits and units for every number and time: "1 hour", "5 GB". Write "may" only for permission. No "soon" or "a
  while" when the code knows the answer.
- R8. One term per thing and one verb per action. The term list is `scripts/text-rules.json`. Write "k-mer", not
  "kmer". A flag or default in the docs matches `--help` word for word.
- R9. No "we" unless the text says who "we" is.
- R10. Every number, limit, default or guarantee matches the code or config that sets it. A change to the value
  changes the text in the same commit.
- R12. A `--help` line fits on 1 line and ends with the default in brackets, such as `[0.7]`.
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

`python3 scripts/check_text_rules.py` checks the files in `scripts/text-rules.json` against
`scripts/text-rules-baseline.json`. It rejects new hits for R1 (semicolons only), R4, R5, R7, R8 and R9. R2 and R6
are advisory in `--report`. Reviewers check em dashes, colons and the other rules. Existing hits are allowed. A file
that loses a hit fails until you lower its baseline in the same commit with `--update-baseline`, which refuses
increases. Mark a deliberate exception with a comment `text-rules-allow: R1 <reason>` on that line or the line
before. Text outside these files, such as commit messages, is not linted. Check it against this list by hand.

