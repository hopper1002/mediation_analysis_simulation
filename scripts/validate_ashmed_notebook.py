"""Check notebook 03's executed R outputs and standalone HTML."""
from pathlib import Path
from html.parser import HTMLParser
import nbformat

root = Path(__file__).resolve().parents[1]
path = root / "notebooks/03_ashmed_comparison.ipynb"
nb = nbformat.read(path, as_version=4)
nbformat.validate(nb)
assert nb.metadata.kernelspec.name == "ir"
codes = [c for c in nb.cells if c.cell_type == "code"]
assert all(c.execution_count is not None for c in codes)
errors = [o for c in codes for o in c.outputs if o.output_type == "error"]
assert not errors, errors
images = sum("image/png" in o.get("data", {}) for c in codes for o in c.outputs)
assert images == 8, images # six four-method main plots and two shrinkage illustrations
source = "\n".join(c.source for c in codes)
assert "generate_data_files(" not in source and "simulate_dataset(" not in source
assert "read_ashmed_inputs(" in source
html = path.with_suffix(".html")
assert html.exists()

class ImageCounter(HTMLParser):
    count = 0
    def handle_starttag(self, tag, attrs):
        if tag == "img" and dict(attrs).get("src", "").startswith("data:image/png"):
            self.count += 1

parser = ImageCounter()
parser.feed(html.read_text(encoding="utf-8"))
assert parser.count == images, (parser.count, images)
print(f"Notebook 03: {len(codes)} R cells executed, no errors, {images} embedded figures; HTML matched; no data generator calls.")
