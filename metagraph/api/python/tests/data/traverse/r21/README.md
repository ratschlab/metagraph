# Level-5 documents for R21 (4)

`merge_level5.graphlet.json` and `annotate_level5.graphlet.json` are the CLI fixtures
`documents/merge.graphlet.json` and `documents/annotate.graphlet.json` as the level-5 walker
wrote them (revision 67bef367, before the owner's decision R21 (4)): their merge at 62
displays the parent that arrived first, which carries fewer labels than the other. The
regenerated documents (feature level 6) display the majority parent there.

`test_traverse_coordinates.TestDisplayedParentAtMerges` compares the two generations of the
same locus: compare() must call a comparison across the rule change 'qualified', never a
definite difference. Nothing regenerates these files; they stay as recorded.
