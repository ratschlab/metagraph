# Engineering communication

Use [the shared principles](PRINCIPLES.md) for meaning, audience and evidence. The formats below serve different decisions. Read only the relevant section. Repository templates, contribution policies and required disclosures take precedence. This guide does not grant permission to publish, approve or merge.

The cited sources are established engineering practice, not experiments establishing an optimal template. Additional choices are labeled **house guidance**. There is no universal word count, required section list or mandatory verdict phrase. These checks require a reader and are not enforced by Vale.

## Pull request descriptions

Help a reviewer understand the problem, resulting behavior and reason for the approach. Write a specific title that remains useful in search. Give enough context for someone who has not seen the discussion; summarize essential reasoning before linking to it. Update the description when the scope changes. **Practice:** Google's change-description guidance emphasizes what changed, why, and the lasting record. [S20](SOURCES.md#engineering-communication)

**House guidance:** lead with the concrete problem and result. For a bug, a small trigger and before/after example often explains more than a file inventory. Include compatibility, operational risk or migration steps when they affect the review. Describe meaningful manual checks and material limitations; avoid duplicating routine CI status or test counts. A small change can need only a paragraph. Include discarded alternatives only when they explain a material tradeoff. Keep agent conversation history and local scratch artifacts out of the description.

Example: “An interrupted upload leaves a partial file under its final name. Write to a temporary file and rename it after completion, so readers only see complete uploads.” Add atomicity or durability claims only after verifying them for the actual platform.

## Reviews and review replies

Explain the problem and why it matters. Discuss the code or behavior respectfully. Distinguish a required correction from a question or optional improvement, and allow the author to choose a suitable fix. **Practice:** Google's review-comment guidance favors reasoned, actionable feedback; its review standard distinguishes code-health concerns from personal preference. [S21, S22](SOURCES.md#engineering-communication)

**House guidance:** attach a finding to the smallest useful location, state the triggering condition and impact, then provide evidence or a concrete example. Label uncertainty as uncertainty. Use the project's severity scale when it has one; do not invent severity to make a preference blocking. Check the current revision and existing discussion before repeating a finding. A suggested fix is useful when clear, but a valid finding need not contain a full redesign.

Example: “If the server commits the request before this timeout, retrying can create a duplicate job. Check the request ID before resubmitting.” This is a hypothesis until supported by the server's retry contract or a reproduction.

For a reply, address the actual concern: explain the correction and point to it, or explain why the current behavior is intentional with evidence. Recheck the revised code before treating the finding as addressed. In a review summary, separate what was inspected, what was demonstrated and any material coverage gap. Do not imply human review, approval or exhaustive coverage that did not occur. Thread resolution and approval authority remain project policy.

## Commit titles and bodies

Make the subject describe the change on its own. In the body, preserve the motivation, prior behavior and important design constraints that a future maintainer cannot readily recover from the diff. Essential reasoning should remain understandable without a transient external discussion. **Practice:** Git's contribution guide provides this model; its detailed formatting requirements are specific to Git. [S23](SOURCES.md#engineering-communication)

**House guidance:** use an imperative subject by default, such as “Reject expired approval requests”. Follow the destination project's prefixes, character limits and trailer rules. Conventional Commits, 50/72-character limits and issue-ID placement are not universal requirements. A body is optional when the subject and diff are sufficient. Describe the individual commit's scope; do not claim the whole feature for a preparatory commit. For squash merges, review the final commit message so it describes the resulting change. Preserve required attribution and disclosures; never invent a sign-off or certification.

## Bug reports and investigation updates

Give the smallest useful reproduction, relevant versions and inputs, expected behavior and actual behavior. Keep observations separate from suspected causes. **Practice:** Simon Tatham's reporting guidance centers on enabling someone else to reproduce or diagnose the failure. [S24](SOURCES.md#engineering-communication)

**House guidance:** keep a short reproducer in the report and link lengthy logs. Remove secrets and unnecessary personal data. State which environment or revision produced the evidence. For intermittent behavior, give the observed frequency and conditions rather than claiming determinism. A successful workaround is not proof that the root cause is fixed.

## Release notes and migration instructions

Explain notable effects for users, grouped by the project's release conventions. **Practice:** Keep a Changelog treats release history as a curated, human-readable record rather than a dump of commit messages. [S25](SOURCES.md#engineering-communication)

**House guidance:** identify affected versions, changed defaults, compatibility breaks and required user action. Where migration is needed, provide an ordered path, prerequisites and a way to verify completion. Describe rollback only when it is supported; disclose irreversible steps. Do not copy internal refactoring details into release notes unless they explain a user-visible effect. Follow the project's security-disclosure process for sensitive fixes.

## Decision requests and handoffs

This section is **house guidance** applying P1, P4 and P7. A decision request names the pending decision, relevant facts, options, recommendation, consequences and exact scope of the requested authorization. Avoid asking again when authorization already covers the action. Do not convert a recommendation or an agent's review into someone else's approval.

A handoff records the objective, current revision, verified outcome, outstanding work and next action. Include reproducible checks and locations of necessary artifacts. For an ongoing job, identify how to check its state. Separate historical observations from current state; do not carry a stale “running” or “approved” claim forward without its scope and evidence. Keep the record only as detailed as the next person needs to continue.

## Agent self-check

Before presenting or sending a draft, read it as the recipient: is the point clear, are the claims supported, and is the next action understandable? Remove repetition and irrelevant process history. Retain uncertainty and limitations that change the decision. Follow the repository's actual approval and disclosure requirements. This is a review aid, not a requirement to append a compliance checklist to every message.
