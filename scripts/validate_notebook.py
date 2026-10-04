"""Validate the delivered notebook structure and outputs, without executing it."""
from pathlib import Path
import nbformat

root = Path(__file__).resolve().parents[1]
path = root / "notebooks" / "01_simulation_workbench.ipynb"
nb = nbformat.read(path, as_version=4)
nbformat.validate(nb)
assert nb.metadata.kernelspec.name == "ir"
code = [c for c in nb.cells if c.cell_type == "code"]
assert all(c.execution_count is not None for c in code)
errors = [o for c in code for o in c.outputs if o.output_type == "error"]
assert not errors, errors
images = [o for c in code for o in c.outputs if "image/png" in o.get("data", {})]
assert len(images) == 3, len(images)
assert path.with_suffix(".html").stat().st_size > 200000
assert not any("Question? Answer." in c.source for c in nb.cells)
streams = "\n".join(o.get("text", "") for c in code for o in c.outputs if o.output_type == "stream")
assert "Analysed 720/720 saved datasets; method errors: 0" in streams
assert "3个固定初始化" in "\n".join(c.source for c in nb.cells)
print(f"PASS: valid R notebook, {len(code)} executed cells, no error outputs, three paired FDR/power figures, HTML export.")
