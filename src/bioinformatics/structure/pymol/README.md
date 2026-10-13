# structure/pymol/

Session-oriented PyMOL recipes (not a Python import package). Open PyMOL and
`run` the script you need.

| File | Typical use |
|------|-------------|
| `imageCartoon.py` | Cartoon rendering / figure |
| `imagePutty.py` | Putty (B-factor / flexibility) view |
| `imagePutty_KDnTail.py` | Putty for kinase domain + tail |
| `orientKinase.py` | Standard kinase orientation |
| `KD_tail.py` | Kinase domain + tail selection helpers |
| `mutateProtein.py` | In-silico mutation helpers |
| `snap.py` | Snapshot / image capture |
| `view_matrix.pse` | Saved PyMOL session / view |

Example (inside PyMOL):

```text
run /home/anup/myScripts/src/bioinformatics/structure/pymol/imageCartoon.py
```

Parent: [`../README.md`](../README.md).
