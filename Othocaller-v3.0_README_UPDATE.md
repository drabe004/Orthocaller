# 🧬 Orthocaller v3

**Orthocaller v3** is a gene tree-based orthogroup classification tool for comparative genomics. It uses reconciled gene trees (e.g., from **GeneRax**) and a rooted species tree to identify conserved, duplicated, and lost orthologs across species.

---

## 🚨 What's New in v3?

This version introduces **flexible in-paralog resolution**, with three user-selectable strategies for deciding which in-paralog copy to retain:

| Strategy            | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `farthest` (default) | Retains the most diverged in-paralog (longest branch from duplication node). |
| `shortest`          | Retains the least diverged in-paralog (shortest branch).                     |
| `average_divergence`| Retains the copy whose branch length ratio is closest to the duplication node — a proxy for balanced divergence. |

---
