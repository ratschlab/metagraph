"""A fake metagraph.traverse: just enough of the library's surface for run_graphlet_local's child process (a store,
budget.WORK_MODEL, parser.standalone_text, mcp_tools.GraphletTools and tool_names), with tools that misbehave on
purpose: they hang, answer what is no JSON, end the process, repeat a cursor, print, read stdin."""


class GraphletStore(object):
    def __init__(self, spool_dir, **kw):
        self.spool_dir = spool_dir
