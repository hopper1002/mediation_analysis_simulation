from pathlib import Path
from html.parser import HTMLParser
import nbformat

root = Path(__file__).resolve().parents[1]
path = root / "notebooks/04_ashmed_adaptive.ipynb"
nb = nbformat.read(path, as_version=4)
nbformat.validate(nb)
assert nb.metadata.kernelspec.name == "ir"
codes = [c for c in nb.cells if c.cell_type == "code"]
assert all(c.execution_count is not None for c in codes)
errors = [o for c in codes for o in c.outputs if o.output_type == "error"]
assert not errors, errors
images = sum("image/png" in o.get("data", {}) for c in codes for o in c.outputs)
assert images == 10, images
source = "\n".join(c.source for c in codes)
assert "generate_data_files(" not in source and "simulate_dataset(" not in source
assert "run_final_adaptive_comparison(" in source
class Counter(HTMLParser):
    count = 0
    def handle_starttag(self, tag, attrs):
        if tag == "img" and dict(attrs).get("src", "").startswith("data:image/png"):
            self.count += 1
counter = Counter()
counter.feed(path.with_suffix(".html").read_text(encoding="utf-8"))
assert counter.count == images
print(f"Notebook 04: {len(codes)} R cells executed; no errors; {images} embedded figures; HTML matched; no data generator calls.")
