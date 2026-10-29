# sample-input

This app consumes the *outputs* of the other immune cytometry steps (fcs-prep, peacoqc, cleanet,
flowsom, cytometry-umap, diffcyt, cytometry-batch-assess, cytonorm, cycombine), so there is no
public file to download for it directly.

For a local run, generate the synthetic fixture (3 batches x 2 conditions x 2 samples, every port):

```sh
make fixture      # writes data/fixture/linked/ and data/fixture/flat/ (gitignored)
make run          # runs the Docker image over data/fixture/linked/
python -m tests.make_fixture data/fixture-a --layout linked --workflow a   # workflow A inputs only
python -m tests.make_fixture data/fixture-b --layout linked --workflow b   # workflow B inputs only
```

For real data, run the upstream steps on the demo dataset (Koenig Z, Andorf S, Tilburgs T,
Eswar S. CytoBatchFlagR: Spectral flow cytometry dataset. Zenodo 2025.
doi:10.5281/zenodo.15388817, CC BY 4.0) and point the report at their outputs, one folder per
port, e.g. with workflow-immune-cytometry's `data/work/`:

```
input/metadata      -> data/work/prep/metadata
input/qc            -> data/work/peacoqc/qc
input/doublets      -> data/work/cleanet/qc
input/clusters      -> data/work/flowsom/clusters
input/embedding     -> data/work/umap
input/da, input/ds  -> data/work/diffcyt/da, data/work/diffcyt/ds
input/assessment/<LABEL> -> data/work/batch-assess/<LABEL>   (raw, cytonorm, cycombine)
input/normalization/<method> -> data/work/<method>/qc/<method>  (cytonorm, cycombine)
```

With Docker, mount each as its own volume (symlinks to host paths do not resolve inside the
container); `INPUT_HOST=<folder> make run` works for a folder of real copies.
