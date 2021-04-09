
SMAKE_BASE="./run_snakemake_example.sh"

run:
	$(SMAKE_BASE)

graph:
	"$(SMAKE_BASE)" --forceall --dag > rule_graph.dot
	dot -Tpng rule_graph.dot > rule_graph.png

shell-cmds:
	"$(SMAKE_BASE)" --forceall --dryrun --printshellcmds | grep '^        [A-Za-z]' | sed -E 's/ +/ /g'
