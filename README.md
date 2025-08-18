# Genesis-Sankoff
A Sankoff-Rousseau-like algorithm for inferring horizontal gene transfers using tree-based networks
This repository contains the implementation of the **Genesis Algorithm** for the inference of **Horizontal Gene Transfers (HGTs)** in phylogenetic networks.  
The work builds upon our [Perfect Transfer Network model](https://github.com/AliLopSan/ptns) and is part of the methods developed for our article: **A Sankoff-Rousseau-like Algorithm for Minimizing Lateral
    Gene Transfers and Losses on Single Origin Characters**.

The project introduces **tree-based networks** as a data structure to capture evolutionary histories where both **vertical inheritance** and **horizontal gene transfer** events occur.  

---

## Key Features

- **Tree-Based Network Class (`TB_Network`)**
  - Represents a phylogenetic tree with additional transfer arcs.
  - Supports efficient traversal methods (preorder, postorder, subtree traversal).
  - Tracks transfer edges and their weights.

- **Extended Sankoff Algorithm**
  - Implements **orign-tracking labeling** with costs for **losses** and **first appearances**.

- **Transfer Inference Tools**
  - Detects **transfer highways** (recurring donor–recipient edges with strong support).
  - Identifies **replaceable clades** and **triangle transfers** for fine-grained evolutionary interpretation.
  - Greedy completion procedure to build transfer scenarios consistent with observed data.

- **Integration with [Tralda](https://github.com/david-schaller/tralda/)**
  - Interfaces with Tralda’s tree structures to initialize base trees and compute LCAs efficiently.

---

## Repository Structure
├── src/ # Core implementation (treebased.py, algorithm code)
├── data/ # Input species trees and gene trees (NHX format)
├── results/ # Inferred transfers, statistics, and network outputs
├── notebooks/ # Jupyter notebooks for reproducible figures & analyses
└── README.md # Project documentation (this file)

