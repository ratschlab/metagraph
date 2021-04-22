
SMAKE_BASE="./run_snakemake_example.sh"

run:
	$(SMAKE_BASE)

graph:
	"$(SMAKE_BASE)" --forceall --dag > example_workflow/rule_graph.dot
	dot -Tpng example_workflow/rule_graph.dot > example_workflow/rule_graph.png

shell-cmds:
	"$(SMAKE_BASE)" --forceall --dryrun --printshellcmds | grep '^        [A-Za-z]' | sed -E 's/ +/ /g'
