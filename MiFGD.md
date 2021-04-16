---
layout: default
---

*(Blog post by Junhyung Lyle Kim)*


This blog post is about my recent work on quantum state tomography using (accelerated) non-convex programming.[^kim2021fast] Our manuscript is on [arXiv](https://arxiv.org/abs/2104.07006). This is a joint work with my advisor [Prof. Tasos Kyrillidis](https://akyrillidis.github.io/about/) at Rice University, [Dr. Amir Kalev](https://scholar.google.com/citations?user=te_1dnAAAAAJ&hl=en) at USC, and [Dr. Georgios Kollias](https://researcher.watson.ibm.com/researcher/view.php?person=us-gkollias) and [Dr. Ken Wei](https://scholar.google.com/citations?user=9uuZX3IAAAAJ&hl=en) at IBM Quantum.

## Introduction

Quantum state tomography (QST) is one of the main procedures to identify the nature of imperfections in quantum processing unit (QPU) implementation. High-level procedure is to measure the quantum system, estimate the density matrix using the measured data, and analyze the "fit" between the estimated density matrix and the true density matrix.

QST is generally not scalable due to two bottlenecks: $$i)$$ large data has to be collected to perform tomography; and $$ii)$$ the space of density matrices grows exponentially, from which the one that is consistent with the data has to be found.

To address the first bottleneck, prior information is often assumed and leveraged to reduce the number of data required. For example, in compressed sensing QST, [^gross2010quantum] [^kalev2015quantum]  it assumes that the density matrix is of low-rank. Similarly, in neural network QST, the wavefunctions are assumed to be real and positive. [^torlai2018neural] [^torlai2019machine] [^beach2019qucumber]

To give a concrete example, in the figure below, real (top) and imaginary (bottom) parts of four different states are shown: $$i)$$ $$\texttt{GHZ}$$ state, $$ii)$$ $$\texttt{GHZminus}$$ state, $$iii)$$ $$\texttt{Hadamard}$$ state, and $$iv)$$ $$\texttt{Random}$$ state; for the mathematical description of the above states, refer to our paper. [^kim2021fast] As can be seen, for $$\texttt{GHZ}$$ and $$\texttt{GHZminus}$$ states, only four corners of the real parts have non-zero entries. Therefore, the density matrices of these states are both of low-rank and sparse. If these kinds of "structures" are smartly leveraged, one can sometimes confine the search space of density matrices greatly, leading to less number of measurements required for successful tomography results. 


[^gross2010quantum]: D. Gross, Y.-K. Liu, S. Flammia, S. Becker, and J. Eisert. Quantum state tomography via compressed
sensing. Physical review letters, 105(15):150401, 2010.

[^kalev2015quantum]: A. Kalev, R. Kosut, and I. Deutsch. Quantum tomography protocols with positivity are compressed
sensing protocols. NPJ Quantum Information, 1:15018, 2015.

[^torlai2018neural]: Giacomo Torlai, Guglielmo Mazzola, Juan Carrasquilla, Matthias Troyer, Roger Melko, and Giuseppe
Carleo. Neural-network quantum state tomography. Nat. Phys., 14:447–450, May 2018.

[^torlai2019machine]: Giacomo Torlai and Roger Melko. Machine-learning quantum states in the NISQ era. Annual Review
of Condensed Matter Physics, 11, 2019.

[^beach2019qucumber]: Matthew JS Beach, Isaac De Vlugt, Anna Golubeva, Patrick Huembeli, Bohdan Kulchytskyy, Xiuzhe
Luo, Roger G Melko, Ejaaz Merali, and Giacomo Torlai. Qucumber: wavefunction reconstruction with neural networks. SciPost Physics, 7(1):009, 2019.

[^goncalves2016projected]: D. Gonçalve, M. Gomes-Ruggiero, and C. Lavor. A projected gradient method for optimization over
density matrices. Optimization Methods and Software, 31(2):328–341, 2016.

[^bolduc2017projected]: E. Bolduc, G. Knee, E. Gauger, and J. Leach. Projected gradient descent algorithms for quantum state tomography. npj Quantum Information, 3(1):44, 2017.

[^shang2017superfast]: Jiangwei Shang, Zhengyun Zhang, and Hui Khoon Ng. Superfast maximum-likelihood reconstruction
for quantum tomography. Phys. Rev. A, 95:062336, Jun 2017.

[^hu2019reconstructing]: Zhilin Hu, Kezhi Li, Shuang Cong, and Yaru Tang. Reconstructing pure 14-qubit quantum states in
three hours using compressive sensing. IFAC-PapersOnLine, 52(11):188 – 193, 2019. 5th IFAC Conference on Intelligent Control and Automation Sciences ICONS 2019.

[^hou2016full]: Zhibo Hou, Han-Sen Zhong, Ye Tian, Daoyi Dong, Bo Qi, Li Li, Yuanlong Wang, Franco Nori, Guo-Yong Xiang, Chuan-Feng Li, et al. Full reconstruction of a 14-qubit state within four hours. New Journal of Physics, 18(8):083036, 2016.

[^kim2021fast]: Junhyung Lyle Kim, George Kollias, Amir Kalev, Ken X. Wei, Anastasios Kyrillidis. Fast quantum state reconstruction via accelerated non-convex programming. arXiv preprint arXiv:2104.07006, 2021.

[^kyrillidis2018provable]: A. Kyrillidis, A. Kalev, D. Park, S. Bhojanapalli, C. Caramanis, and S. Sanghavi. Provable quantum state tomography via non-convex methods. npj Quantum Information, 4(36), 2018.

[^recht2010guaranteed]: Benjamin Recht, Maryam Fazel, and Pablo A Parrilo. Guaranteed minimum-rank solutions of linear matrix equations via nuclear norm minimization. SIAM review, 52(3):471–501, 2010.




[back](./)
