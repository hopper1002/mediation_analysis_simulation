"""Check the delivered teaching notebook, embedded figures and HTML export."""
from pathlib import Path
import base64
from html.parser import HTMLParser
import nbformat

root = Path(__file__).resolve().parents[1]
path = root/"notebooks/02_teaching_mediation.ipynb"
nb = nbformat.read(path,as_version=4)
nbformat.validate(nb)
assert nb.metadata.kernelspec.name == "ir"
code = [c for c in nb.cells if c.cell_type=="code"]
assert all(c.execution_count is not None for c in code)
assert not [o for c in code for o in c.outputs if o.output_type=="error"]
images = [o.data["image/png"] for c in code for o in c.outputs if "image/png" in o.get("data",{})]
assert len(images)==5, len(images)
assert all(base64.b64decode(x).startswith(b"\x89PNG\r\n\x1a\n") for x in images)
streams = "\n".join(o.get("text","") for c in code for o in c.outputs if o.output_type=="stream")
assert "Analysed 720/720 saved datasets; method errors: 0" in streams
assert "720 / 720" in streams
sources = "\n".join(c.source for c in nb.cells)
assert sources.index("generated <- generate_data_files") < sources.index("experiment <- run_experiment")
assert "current_teaching.rds" in sources and "config$repetitions <- 40L" in sources
html = path.with_suffix(".html").read_text(encoding="utf-8")
class ImageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []
    def handle_starttag(self, tag, attrs):
        src = dict(attrs).get("src", "")
        if tag == "img" and src.startswith("data:image/png;base64,"):
            self.images.append(base64.b64decode(src.split(",",1)[1]))
parser = ImageParser()
parser.feed(html)
# nbconvert also embeds a small PNG inside CSS; count actual figure img tags.
assert len(html)>200000 and parser.images == [base64.b64decode(x) for x in images]
print(f"PASS: valid R teaching notebook; {len(code)} executed cells; no errors; five embedded PNGs; standalone HTML; data generated before analysis.")
