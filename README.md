# Autoseed v2.2 (notebook + pip version)

This package installs:

- A Jupyter notebook interface (`autoseed_v2.2.ipynb`)
- Four compiled C programs:
  - `localmax-motif`
  - `seedextender`
  - `motifsimilarity`
  - `genint-PWM`

During installation (`pip install .` or via GitHub), the C source is downloaded
from the original repositories and compiled automatically.

## Installation (from GitHub)

```bash
pip install git+https://github.com/jutaipal/autoseed-notebook.git
