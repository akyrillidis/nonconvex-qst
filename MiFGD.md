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


{{< figure src="./assets/img/state-plots.png" title="From left to right: GHZ, GHZminus, Hadamard, and Random states. All states are in 4-qubit system. " width="100%">}}


With regards to the second bottleneck, variants of gradient descent convex solvers were proposed under synthetic scenarios. [^goncalves2016projected] [^bolduc2017projected] [^shang2017superfast] [^hu2019reconstructing] However, due to the exponentially increasing space of density matrices, these methods often can be only applied to relatively small system, on top of relying on special-purpose hardwares and proper distributed system designs.[^hou2016full]  

On the contrary, non-convex optimization methods can perform much faster. It was recently shown that one can formulate compressed sensing QST as a non-convex problem,[^kyrillidis2018provable] which can be solved with rigorous convergence guarantees, allowing density matrix estimation in a large system. A relevant result can be seen in the **Results** section below, where we compare our proposed (accelerated) non-convex method with [convex methods from $$\texttt{Qiskit}$$](https://qiskit.org/documentation/stubs/qiskit.ignis.verification.TomographyFitter.html).  

In this work, we consider the setup where $$n$$-qubit state is close to a pure state, thus its density matrix is of low-rank. We introduce an accelerated non-convex algorithm with provable gaurantees, which we call $$\texttt{MiFGD}$$, short for "$$\texttt{M}$$omentum $$\texttt{i}$$nspired $$\texttt{F}$$actored $$\texttt{G}$$radient $$\texttt{D}$$escent."

## Problem setup
We consider the reconstruction of a low-rank density matrix $$\rho^\star \in \mathbb{C}^{d \times d}$$ on a $$n$$-qubit Hilbert space, where $$d=2^n$$, through the following $$\ell_2$$-norm reconstruction objective:

\begin{align}
\label{eq:objective} \tag{1}
\min_{\rho \in \mathbb{C}^{d \times d}}
\quad & f(\rho) := \tfrac{1}{2} ||\mathcal{A}(\rho) - y||_2^2 \quad
\text{subject to} \quad \rho \succeq 0, ~\texttt{rank}(\rho) \leq r.
\end{align}

Here, $$y \in \mathbb{R}^m$$ is the measured data through quantum computer or simulation, and $$\mathcal{A}(\cdot): \mathbb{C}^{d \times d} \rightarrow \mathbb{R}^m$$ is the linear sensing map. The sensing map relates the density matrix $$\rho$$ to the measurements through [Born rule](https://en.wikipedia.org/wiki/Born_rule): $$\left( \mathcal{A}(\rho) \right)_i = \text{Tr}(A_i \rho),$$ where $$A_i \in \mathbb{C}^{d \times d},~i=1, \dots, m$$ are the sensing matrices. From the objective function above, we see two constraints: $$i)$$ the density matrix $$\rho$$ is a positive semidefinite matrix (which is a convex constraint), and $$ii)$$ the rank of the density matrix is less than $$r$$ (which is a non-convex constraint). 

As mentioned earlier, we focus on _compressed sensing quantum state tomography_ setting, where the number of measured data $$m$$ is much smaller than the problem dimension $$O(d^2)$$. Compressed sensing is a powerful optimization framework developed mainly by [Emmanuel Candès](https://statweb.stanford.edu/~candes/), [Justin Romberg](https://jrom.ece.gatech.edu/), [Terence Tao](https://www.math.ucla.edu/~tao/) and [David Donoho](https://web.stanford.edu/dept/statistics/cgi-bin/donoho/), and requires the following pivotal assumption on the sensing matrix $$\mathcal{A}(\cdot)$$, namely the **Restricted Isometry Property (RIP)** (on $$\texttt{rank}$$-$$r$$ matrices): [^recht2010guaranteed]

\begin{align}
\label{eq:rip} \tag{2}
(1 - \delta_r) \cdot  || X ||_F^2 \leq || \mathcal{A}(X) ||_2^2 \leq (1 + \delta_r) \cdot ||X||_F^2.
\end{align}

Intuitively, the above RIP assumption states that the sensing matrices $$\mathcal{A}(\cdot)$$ only "marginally" changes the norm of the matrix $$X$$.

Going back to the main optimization problem in Eq. \eqref{eq:objective}, instead of solving it as is, we propose to solve a factorized version of it, following recent work [^kyrillidis2018provable]:
\begin{align}
\label{eq:factored-obj} \tag{3}
  \min_{U \in \mathbb{C}^{d \times r}} f(UU^\dagger) := \tfrac{1}{2} || \mathcal{A} (UU^\dagger) - y ||_2^2,
\end{align}
where $$U^\dagger$$ denotes the [adjoint](https://en.wikipedia.org/wiki/Conjugate_transpose) of $$U$$. The motivation is rather clear: in the original objective function in Eq. \eqref{eq:objective}, the density matrix $$\rho$$ is represented as a $$d \times d$$ Hermitian matrix, and due to the (non-convex) $$\texttt{rank}(\cdot)$$ constraint, some method to project onto the set of low-rank matrices is required. Instead, we work in the space of the factors $$U \in \mathbb{C}^{d \times r}$$, and by taking an outer-product, the $$\texttt{rank}(\cdot)$$ constraint and the PSD constraint $$\rho \succeq 0$$ are automatically satisfied, leading to the non-convex formulation in Eq. \eqref{eq:factored-obj}. But how do we solve such problem?

A common approach is to use gradient descent on $$U$$, which iterates as follows:
\begin{align}
\label{eq:fgd} \tag{4}
U_{i+1} &= U_{i} - \eta \nabla f(U_i U_i^\dagger) \cdot U_i \\\\
&= U_{i} - \eta \mathcal{A}^\dagger \left(\mathcal{A}(U_i U_i^\dagger) - y\right) \cdot U_i.
\end{align}
Here, $$\mathcal{A}^\dagger: \mathbb{R}^m \rightarrow \mathbb{C}^{d \times d}$$ is the adjoint of $$\mathcal{A}$$, defined as $$\mathcal{A}^\dagger = \sum_{i=1}^m x_i A_i.$$ $$\eta$$ is a hyperparameter called step size or learning rate. This method is called "$$\texttt{F}$$actored $$\texttt{G}$$radient $$\texttt{D}$$escent" ($$\texttt{FGD}$$), and was utilized to solve the non-convex objective function in Eq. \eqref{eq:factored-obj}, (surprisingly) with provable gaurantees.[^kyrillidis2018provable]

## Momentum-inspired Factored Gradient Descent
Momentum is one of the de facto techniques to achieve acceleration in gradient descent type of algorithms. Acceleration methods exist in various forms, including Polyak's momentum, Nesterov's acceleration, classical momentum, etc. They end up behaving pretty similarly, and we will not get into the detail of different acceleration methods in this post. For interested readers, I recommend this [blog post](https://jlmelville.github.io/mize/nesterov.html) by James Melville.

A common feature accross acceleration methods is that, with proper hyper-parameter tuning, they can provide accelerated convergence rate with virtually no additional computation. This is exactly the motivation of the $$\texttt{MiFGD}$$ algorithm we propose for solving the non-convex objective in Eq. \eqref{eq:factored-obj}, and the algorithm proceeds as follows:
\begin{align}
\label{eq:mifgd} \tag{5}
U_{i+1} &= Z_{i} - \eta \mathcal{A}^\dagger \left(\mathcal{A}(Z_i Z_i^\dagger) - y\right) \cdot Z_i, \quad 
Z_{i+1} &= U_{i+1} + \mu \left(U_{i+1} - U_i\right). 
\end{align}

Here, $$Z_i \in \mathbb{C}^{d\times r}$$ is a rectangular matrix (with the same dimension as $$U_i$$) which accumulates the "momentum" of the iterates $$U_i$$. $$\mu$$ is the momentum parameter that balances the weight between the previous estimate $$U_i$$ and the current estimate $$U_{i+1}.$$

While the $$\texttt{MiFGD}$$ algorithm in Eq. \eqref{eq:mifgd} looks quite similar to $$\texttt{FGD}$$ in Eq. \eqref{eq:fgd}, it complicates the convergence theory significantly. This is because the two-step momentum procedure has to be considered, on top of the fact that the objective function in Eq. \eqref{eq:factored-obj} is non-convex. We will not get into the details of the convergence thoery here; interested readers are referred to our paper.[^kim2021fast] We finish this section with an informal theorem that illustrates the convergence behavior of $$\texttt{MiFGD}$$:

**Theorem 1** ($$\texttt{MiFGD}$$ convergence rate (informal)). Assume that $$\mathcal{A}(\cdot)$$ satisfies the RIP for some constant $$0 < \delta_{2r} <1$$. Let $$y = \mathcal{A}(\rho^\star)$$ denote the set of measurements, by measuring the quantum density matrix $$\rho^\star$$. Given a good initialization point $$U_0$$, and setting step size $$\eta$$ and momentum $$\mu$$ appropriately, $$\texttt{MiFGD}$$ converges with a linear rate to a region—with radius that depends on $$O(\mu)$$—around the global solution $$\rho^\star$$. 




[^gross2010quantum]: D. Gross, Y.-K. Liu, S. Flammia, S. Becker, and J. Eisert. Quantum state tomography via compressed sensing. Physical review letters, 105(15):150401, 2010.

[^kalev2015quantum]: A. Kalev, R. Kosut, and I. Deutsch. Quantum tomography protocols with positivity are compressed sensing protocols. NPJ Quantum Information, 1:15018, 2015.

[^torlai2018neural]: Giacomo Torlai, Guglielmo Mazzola, Juan Carrasquilla, Matthias Troyer, Roger Melko, and Giuseppe Carleo. Neural-network quantum state tomography. Nat. Phys., 14:447–450, May 2018.

[^torlai2019machine]: Giacomo Torlai and Roger Melko. Machine-learning quantum states in the NISQ era. Annual Review of Condensed Matter Physics, 11, 2019.

[^beach2019qucumber]: Matthew JS Beach, Isaac De Vlugt, Anna Golubeva, Patrick Huembeli, Bohdan Kulchytskyy, Xiuzhe Luo, Roger G Melko, Ejaaz Merali, and Giacomo Torlai. Qucumber: wavefunction reconstruction with neural networks. SciPost Physics, 7(1):009, 2019.

[^goncalves2016projected]: D. Gonçalve, M. Gomes-Ruggiero, and C. Lavor. A projected gradient method for optimization over density matrices. Optimization Methods and Software, 31(2):328–341, 2016.

[^bolduc2017projected]: E. Bolduc, G. Knee, E. Gauger, and J. Leach. Projected gradient descent algorithms for quantum state tomography. npj Quantum Information, 3(1):44, 2017.

[^shang2017superfast]: Jiangwei Shang, Zhengyun Zhang, and Hui Khoon Ng. Superfast maximum-likelihood reconstruction for quantum tomography. Phys. Rev. A, 95:062336, Jun 2017.

[^hu2019reconstructing]: Zhilin Hu, Kezhi Li, Shuang Cong, and Yaru Tang. Reconstructing pure 14-qubit quantum states in three hours using compressive sensing. IFAC-PapersOnLine, 52(11):188 – 193, 2019. 5th IFAC Conference on Intelligent Control and Automation Sciences ICONS 2019.

[^hou2016full]: Zhibo Hou, Han-Sen Zhong, Ye Tian, Daoyi Dong, Bo Qi, Li Li, Yuanlong Wang, Franco Nori, Guo-Yong Xiang, Chuan-Feng Li, et al. Full reconstruction of a 14-qubit state within four hours. New Journal of Physics, 18(8):083036, 2016.

[^kim2021fast]: Junhyung Lyle Kim, George Kollias, Amir Kalev, Ken X. Wei, Anastasios Kyrillidis. Fast quantum state reconstruction via accelerated non-convex programming. arXiv preprint arXiv:2104.07006, 2021.

[^kyrillidis2018provable]: A. Kyrillidis, A. Kalev, D. Park, S. Bhojanapalli, C. Caramanis, and S. Sanghavi. Provable quantum state tomography via non-convex methods. npj Quantum Information, 4(36), 2018.

[^recht2010guaranteed]: Benjamin Recht, Maryam Fazel, and Pablo A Parrilo. Guaranteed minimum-rank solutions of linear matrix equations via nuclear norm minimization. SIAM review, 52(3):471–501, 2010.




[back](./)
