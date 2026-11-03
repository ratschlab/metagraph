# Sources and limits

Sources checked on 2026-11-03. This is a selected primary-source review, not a systematic literature review. Research reports, standards and practitioner guidance answer different questions. No cited study proves that this particular package improves agent output. Text here is original synthesis; linked works retain their own licenses.

## Research

- **S18. John M. Carroll, Penny L. Smith-Kerker, James R. Ford and Sandra A. Mazur-Rimetz (1987), [The Minimal Manual](https://doi.org/10.1207/s15327051hci0302_2), Human–Computer Interaction 3(2), 123–153.** [Full paper](https://swcarpentry.github.io/swc-releases/2017.02/instructor-training/files/papers/carroll-minimal-manual-1987.pdf) reports two experiments on self-instruction for a word processor. Supports task-oriented learning, guided exploration and error recovery in that setting; does not establish that shorter is always better or that reference material should be omitted.
- **S19. Harald Weinreich, Hartmut Obendorf, Eelco Herder and Matthias Mayer (2008), [Not Quite the Average: An Empirical Study of Web Use](https://www.eelcoherder.com/images/publications/2008/not_quite_the_average.pdf), ACM Transactions on the Web 2(1), article 5, DOI 10.1145/1326561.1326566.** A 25-participant longitudinal browsing study documents rapid browsing and differences between users and site types. Supports designing for scanning, with important limits: observational behavior is not a universal reading percentage, causal proof of a word limit, or a measurement of today's interfaces.

- **S1. John Morkes and Jakob Nielsen (1997), [Concise, Scannable, and Objective: How to Write for the Web](https://www.nngroup.com/articles/concise-scannable-and-objective-how-to-write-for-the-web/).** Reports experiments comparing web-writing variants. Supports concise, scannable, objective presentation in that setting. Its measured gains are not transferable predictions for a CLI, scientific proof or agent-generated text. Does not establish our old numeric sentence/paragraph limits.
- **S2. Hoa Loranger and Kate Moran (2017), [Writing Digital Copy for Domain Experts](https://www.nngroup.com/articles/writing-domain-experts/).** Reports usability work with specialist readers. Supports scanability alongside factual depth and verifiability. A practitioner research report, not evidence that all experts want minimal detail.

## Documentation and language

- **S3. Daniele Procida, [Diátaxis](https://diataxis.fr/).** An established documentation framework separating tutorials, how-to guides, reference and explanation. Supports organizing by need; does not specify punctuation or mandatory page lengths.
- **S4. Google, [Developer documentation style highlights](https://developers.google.com/style/highlights).** Maintained technical-writing practice covering audience, tone, language and structure. A style guide, not an experimental result.
- **S5. Microsoft, [Writing tips for global communications](https://learn.microsoft.com/en-us/style-guide/global-communications/writing-tips).** Supports simple syntax, consistency and localization. Microsoft's [semicolon guidance](https://learn.microsoft.com/en-us/style-guide/punctuation/semicolons) prefers simplification but explicitly allows some uses. It does not support a universal punctuation ban.

## Errors and accessibility

The W3C links below are explanatory Understanding documents, which identify and link the normative success criteria. Applying a few writing rules does not establish full WCAG conformance.

- **S6. W3C, [WCAG 2.2 SC 3.3.3 Error Suggestion, AA](https://www.w3.org/WAI/WCAG22/Understanding/error-suggestion.html).** Known correction suggestions for detected input errors, unless security or purpose would be jeopardized. Scope is input assistance, not every possible runtime failure.
- **S7. GOV.UK Design System, [Error message](https://design-system.service.gov.uk/components/error-message/).** Government-service design practice for validation errors. Supports specific correction guidance; its stylistic conventions need context when used elsewhere.
- **S8. Tim Neusesser and Evan Sunwall (2023), [Error-Message Guidelines](https://www.nngroup.com/articles/error-message-guidelines/).** Practitioner guidance on visible, constructive errors that respect user effort. Does not require invented causes or every diagnostic code to be hidden.
- **S9. W3C, [SC 2.4.6 Headings and Labels, AA](https://www.w3.org/WAI/WCAG22/Understanding/headings-and-labels.html).** Describe topic or purpose; no blanket ban on questions.
- **S10. W3C, [SC 2.4.4 Link Purpose (In Context), A](https://www.w3.org/WAI/WCAG22/Understanding/link-purpose-in-context.html).** Context can establish a link's purpose. Standalone descriptive labels are a useful preference, stricter than the minimum criterion.

## Interface presentation

- **S11. Jakob Nielsen (2006), [Progressive Disclosure](https://www.nngroup.com/articles/progressive-disclosure/).** Established UX guidance for initially exposing important functions and deferring specialized ones. Its advice is not a mandate to hide all technical information or force a single fold everywhere.
- **S12. GOV.UK Design System, [Details](https://design-system.service.gov.uk/components/details/).** Avoid hiding information most users need. The [notification-banner guidance](https://design-system.service.gov.uk/components/notification-banner/) discourages excessive banners; this is a design-system convention, not a universal count constraint.
- **S13. IBM Carbon, [Empty states](https://www.carbondesignsystem.com/building-blocks/core/patterns/empty-states).** Distinguishes empty-state situations and constructive next steps. Does not justify banning every table cell that says “none”.
- **S14. W3C, [SC 4.1.3 Status Messages, AA](https://www.w3.org/WAI/WCAG22/Understanding/status-messages.html).** Qualifying status changes must be programmatically determinable without taking focus; not every DOM update should be announced.
- **S15. W3C, [SC 1.4.10 Reflow, AA](https://www.w3.org/WAI/WCAG22/Understanding/reflow.html).** Includes exceptions for content requiring two-dimensional layout. Do not rewrite this as an unconditional ban on horizontal scrolling.
- **S16. W3C, [SC 1.4.12 Text Spacing, AA](https://www.w3.org/WAI/WCAG22/Understanding/text-spacing.html).** Requires tolerating specified user spacing changes without loss of content/function. It does not mandate those values as the site's default styles.
- **S17. W3C, [SC 1.4.8 Visual Presentation, AAA](https://www.w3.org/WAI/WCAG22/Understanding/visual-presentation.html).** Covers adjustable text presentation, including line length. It is not an AA requirement or a requirement to use a particular CSS unit.

## Existing implementations

[Vale's package catalog](https://vale.sh/explorer) includes Google, Microsoft, write-good, proselint, RedHat and alex. Reuse selected upstream rules where their semantics match the project. Do not enable multiple complete styles blindly: their terminology and editorial conventions can conflict. The package's current YAML retains the earlier application's rule IDs for migration; that is compatibility work, not a claim to replace those guides.

This review recommends Google as the first technical-documentation style to trial, Microsoft for products already following its conventions, and selective proselint/write-good checks as suggestions. A whole-corpus comparison is still required before substituting those implementations. Semantic truth, task fit and rendered UI behavior remain outside regex linting.
