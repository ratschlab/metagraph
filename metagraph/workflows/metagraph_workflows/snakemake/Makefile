
SMAKE_BASE="./run_snakemake_test.sh"

run:
	$(SMAKE_BASE)

graph:
	"$(SMAKE_BASE)" --forceall --dag > test_workflow/rule_graph.dot
	dot -Tpng test_workflow/rule_graph.dot > test_workflow/rule_graph.png

shell-cmds:
	"$(SMAKE_BASE)" --forceall --dryrun --printshellcmds | grep '^        [A-Za-z]' | sed -E 's/ +/ /g'
