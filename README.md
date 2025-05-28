# 📌 comp‑bio‑pairwise‑sequence‑alignment

A lightweight C++17 implementation of classic **pairwise sequence alignment** algorithms, packaged for quick experimentation on curated DNA datasets.

---

## 🚀 Features

* **Optimised**: `-O2` build flags and cache‑friendly data structures for competitive runtimes.
* **Modern C++17**: clean, header‑only core – no external dependencies.
* **Reproducible**: datasets bundled under `datasets/`; deterministic output in `output.txt`.

---

## 🛠️ Requirements

| Tool                | Version                          |
| ------------------- | -------------------------------- |
| `g++`               | **9.0** or newer (C++17 support) |
| `make` *(optional)* | any                              |

> 💡 *Clang works too – swap `g++` for `clang++` in the commands below.*

---

## ⚙️ Build & Run

```bash
# compile
$ g++ -O2 -std=c++17 q1_solution.cpp -o q1_solution

# execute – outputs alignment metrics and writes results to output.txt
$ ./q1_solution
```

A one‑liner:

```bash
$ g++ -O2 -std=c++17 q1_solution.cpp -o q1_solution && ./q1_solution
```

---

## 📂 Project Layout

```text
.
├── datasets/        # input FASTA / text sequences
│   ├── NM_013096.txt
│   └── ...
├── output.txt       # program output
├── q1_solution.cpp  # main source file
└── README.md        # you are here
```

---

## 📈 Dataset

The `datasets/` directory ships with curated sequences referenced in the homework hand‑out.

---

## 📝 Output

Alignment statistics are printed to **stdout**.

```text
seq1 (./datasets/NM_013096.txt) length: 557
seq2 (./datasets/NM_008218.txt) length: 569
seq3 (./datasets/NM_000558.txt) length: 577
Total running time (seconds): 3.52415
Peak RSS ≈ 8 MB
Alignment length of the sequences : 592
# of exact matches  (identical XYZ) : 425
Total score : 5284
...
```

---


## 📚 References

* D. Gusfield – *Algorithms on Strings, Trees, and Sequences* (1997)
* S. B. Needleman & C. D. Wunsch – “A general method applicable to the search for similarities in the amino acid sequence of two proteins” (1970)
* T. F. Smith & M. S. Waterman – “Identification of common molecular subsequences” (1981)

---
