"""Execute the R notebook and export a standalone HTML using nbclient/nbconvert."""
import os
import asyncio
from pathlib import Path
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

root = Path(__file__).resolve().parents[1]
if os.name == "nt":
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
path = root / "notebooks" / "01_simulation_workbench.ipynb"
# The host has Unix locale strings that R on Windows rejects during startup.
# This changes only this process and its R-kernel child environment.
for name in ("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG"):
    if os.environ.get(name, "").endswith(".UTF-8"):
        os.environ.pop(name, None)
nb = nbformat.read(path, as_version=4)
client = NotebookClient(nb, timeout=7200, kernel_name="ir", resources={"metadata": {"path": str(root)}})
try:
    client.execute()
finally:
    nbformat.write(nb, path)
body, _ = HTMLExporter().from_notebook_node(nb)
body = body.replace("</head>", """<style>
.jp-RenderedMarkdown table { width:100%; table-layout:fixed; }
.jp-RenderedMarkdown th:first-child { width:18%; }
.jp-RenderedMarkdown td { overflow-wrap:anywhere; }
.jp-RenderedHTMLCommon table { font-size:13px; }
</style></head>""")
html = path.with_suffix(".html")
html.write_text(body, encoding="utf-8")
print(f"Executed {path}", flush=True)
print(f"HTML {html}", flush=True)
print(f"Code cells: {sum(c.cell_type == 'code' for c in nb.cells)}; error outputs: "
      f"{sum(o.output_type == 'error' for c in nb.cells if c.cell_type == 'code' for o in c.outputs)}", flush=True)
