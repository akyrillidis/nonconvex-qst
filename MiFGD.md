---
layout: default
---


This blog post is about my recent work on quantum state tomography using (accelerated) non-convex programming.[^kim2021fast] Our manuscript is on [arXiv](https://arxiv.org/abs/2104.07006). This is a joint work with my advisor [Prof. Tasos Kyrillidis](https://akyrillidis.github.io/about/) at Rice University, [Dr. Amir Kalev](https://scholar.google.com/citations?user=te_1dnAAAAAJ&hl=en) at USC, and [Dr. Georgios Kollias](https://researcher.watson.ibm.com/researcher/view.php?person=us-gkollias) and [Dr. Ken Wei](https://scholar.google.com/citations?user=9uuZX3IAAAAJ&hl=en) at IBM Quantum.

## Introduction

Quantum state tomography (QST) is one of the main procedures to identify the nature of imperfections in quantum processing unit (QPU) implementation. High-level procedure is to measure the quantum system, estimate the density matrix using the measured data, and analyze the "fit" between the estimated density matrix and the true density matrix.

QST is generally not scalable due to two bottlenecks: $$i)$$ large data has to be collected to perform tomography; and $$ii)$$ the space of density matrices grows exponentially, from which the one that is consistent with the data has to be found.

To address the first bottleneck, prior information is often assumed and leveraged to reduce the number of data required. For example, in compressed sensing QST, [^gross2010quantum] [^kalev2015quantum]  it assumes that the density matrix is of low-rank. Similarly, in neural network QST, the wavefunctions are assumed to be real and positive. [^torlai2018neural] [^torlai2019machine] [^beach2019qucumber]


[^kim2021fast]: Junhyung Lyle Kim, George Kollias, Amir Kalev, Ken X. Wei, Anastasios Kyrillidis. Fast quantum state reconstruction via accelerated non-convex programming. arXiv preprint arXiv:2104.07006, 2021.



[back](./)
