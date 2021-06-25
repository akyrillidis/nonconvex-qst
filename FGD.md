---
layout: default
---

This blog post is about our previous work on quantum state tomography using non-convex programming.[^kyrillidis2018provable] This manuscript is on [arXiv](https://arxiv.org/pdf/1711.02524.pdf), but also published at [npj Quantum Information](https://www.nature.com/articles/s41534-018-0080-4) This is a joint work of [Prof. Tasos Kyrillidis](https://akyrillidis.github.io/about/) at Rice University, [Dr. Amir Kalev](https://scholar.google.com/citations?user=te_1dnAAAAAJ&hl=en) at USC, [Dr. Dohyung Park](), [Dr. Srinadh Bhojanapalli](https://bsrinadh.github.io/), [Dr. Constantine Caramanis](https://caramanis.github.io/) and [Dr. Sujay Sanghavi](http://users.ece.utexas.edu/~sanghavi/).

## Introduction

Like any other processor, the behavior of a quantum information processor must be characterized, verified, and certified. Quantum state tomography (QST) is one of
the main tools for that purpose [1]. Yet, it is generally an inefficient procedure, since the number of parameters that specify quantum states, grows exponentially with the number of sub-systems. This inefficiency has two practical manifestations: (i) without any prior information, a vast number of data points needs to be collected [1]; (ii) once the data is gathered, a numerical procedure should be executed on an exponentially-high dimensional space, in order to infer the quantum state that is most consistent with the observations. Thus, to perform QST on nowadays steadily growing quantum processors [2, 3], we must introduce novel, more efficient, techniques for its completion. QST is generally not scalable due to two bottlenecks: $$i)$$ large data has to be collected to perform tomography; and $$ii)$$ the space of density matrices grows exponentially, from which the one that is consistent with the data has to be found.

To improve the efficiency of QST, we need to complement it with numerical algorithms that can efficiently handle large search spaces using limited amount of data, while having rigorous performance guarantees. This is the purpose of this work. Inspired by the recent advances on finding the global minimum in non-convex problems [24–38], we propose the application of alternating gradient descent in QST, that operates directly on the assumed low-rank structure of the density matrix. The algorithm –named Projected Factored Gradient Decent (ProjFGD) and described below in detail– is based on the recently analyzed non-convex method in [29] for PSD matrix factorization problems. The added twist is the inclusion of further constraints in the optimization program, that makes it applicable for tasks such as QST.

## Problem setup
We begin by describing  the problem of QST. 
We are focusing here on  QST of a low-rank $n$-qubit state, $\rho_{\star}$, from measuring expectation values of $n$-qubit Pauli observables $\{P_i\}_{i=1}^m$. We denote by $y \in \mathbb{R}^m$  the measurement vector with elements $y_i = \tfrac{2^n}{\sqrt{m}}\text{Tr}(P_i \cdot \rho_\star)+e_i,~i = 1, \dots, m$, for some measurement error $e_i$. 
The normalization $\tfrac{2^n}{\sqrt{m}}$ is chosen to follow the results of Liu\cite{liu2011universal}.
For brevity, we denote $\linmap : \mathbb{C}^{2^n \times 2^n} \rightarrow \mathbb{R}^m$ as the linear ``sensing" map, such that $(\linmap(\rho))_i = \tfrac{2^n}{\sqrt{m}} \text{Tr}(P_i \cdot \rho)$, for $i = 1, \dots, m$. 

An $n$-qubit Pauli observable is given by $P=\otimes_{j=1}^n s_j$ where $s_j\in\{\mathbb{1},\sigma_x,\sigma_y,\sigma_z\}$.
There are $4^n$ such observables in total. 
In general, one needs to have the expectation values of all $4^n$ Pauli observables to uniquely reconstruct $\rho_\star$. However, since according to our assumption $\rho_\star$ is a low-rank quantum state, we can apply the CS result~\citep{gross2010quantum, liu2011universal},  that guarantees a robust estimation, with high probability, from the measurement of the expectation values of just $m={\cal O}(r 2^n n^6)$ randomly chosen Pauli observables, where $r\ll 2^n$ is the rank of $\rho_\star$.


As mentioned earlier, we focus on _compressed sensing quantum state tomography_ setting, where the number of measured data $$m$$ is much smaller than the problem dimension $$O(d^2)$$. Compressed sensing is a powerful optimization framework developed mainly by [Emmanuel Candès](https://statweb.stanford.edu/~candes/), [Justin Romberg](https://jrom.ece.gatech.edu/), [Terence Tao](https://www.math.ucla.edu/~tao/) and [David Donoho](https://web.stanford.edu/dept/statistics/cgi-bin/donoho/), and requires the following pivotal assumption on the sensing matrix $$\mathcal{A}(\cdot)$$, namely the **Restricted Isometry Property (RIP)** (on $$\texttt{rank}$$-$$r$$ matrices): [^recht2010guaranteed]

\begin{align}
\label{eq:rip} \tag{2}
(1 - \delta_r) \cdot  || X ||_F^2 \leq || \mathcal{A}(X) ||_2^2 \leq (1 + \delta_r) \cdot ||X||_F^2.
\end{align}

Intuitively, the above RIP assumption states that the sensing matrices $$\mathcal{A}(\cdot)$$ only "marginally" changes the norm of the matrix $$X$$.

An accurate estimation of $\rhoo$ is obtained  by solving, essentially, a convex optimization problem constrained to the set of quantum states~\citep{kalev2015quantum}, consistent with the measured data.  
Among the various problem formulations for QST, 
two convex program examples are the trace-minimization program that is typically studied in the context of CS QST:
\begin{equation}
	\begin{aligned}
		& \underset{\rho \in \mathbb{C}^{2^n \times 2^n}}{\text{minimize}}
		& & \text{Tr}(\rho) \\
		& \text{subject to}
		& & \rho \succeq 0, \\
		& & & \|y - \mathcal{M}(\rho)\|_2 \leq \epsilon,
	\end{aligned} \label{eq:CVX1}
\end{equation}
and the least-squares program,
\begin{equation}
	\begin{aligned}
		& \underset{\rho \in \mathbb{C}^{2^n \times 2^n}}{\text{minimize}}
		& & \tfrac{1}{2} \cdot \|y - \mathcal{M}(\rho)\|_2^2 \\
		& \text{subject to}
		& & \rho \succeq 0, \\
		& & & \text{Tr}(\rho) \leq 1,
	\end{aligned} \label{eq:CVX2}
\end{equation}
which is closely related to the (negative) log-likelihood minimization under Gaussian noise assumption.
The constraint $\rho \succeq 0$ captures the positive semi-definite assumption, $\|\cdot\|_2$ is the vector Euclidean $\ell_2$-norm, and $\epsilon >0$ is a parameter related to the error level in the model.
Key in both programs is the combination of the PSD constraint and the trace object: combined, they constitute the tightest convex relaxation to the low-rank, PSD structure of the unknown $\rho_\star$; see also Recht \emph{et al.}\citep{recht2010guaranteed}.  
The constraint $\text{Tr}(\rho) = 1$ is relaxed in~\eqref{eq:CVX2} to allow more robustness to noise, following Kalev \emph{et al.}\citep{kalev2015quantum}. The solutions of these programs should be normalized to have unit trace to represent quantum states.
We note that if $\mathcal{M}$ corresponds to a positive-operator valued measure (POVM), or includes the identity operator, then the explicit trace constraint is redundant.

## Projected Factored Gradient Descent
Momentum is one of the de facto techniques to achieve acceleration in gradient descent type of algorithms. Acceleration methods exist in various forms, including Polyak's momentum, Nesterov's acceleration, classical momentum, etc. They end up behaving pretty similarly, and we will not get into the detail of different acceleration methods in this post. For interested readers, I recommend this [blog post](https://jlmelville.github.io/mize/nesterov.html) by James Melville.

A common feature accross acceleration methods is that, with proper hyper-parameter tuning, they can provide accelerated convergence rate with virtually no additional computation. This is exactly the motivation of the $$\texttt{MiFGD}$$ algorithm we propose for solving the non-convex objective in Eq. \eqref{eq:factored-obj}, and the algorithm proceeds as follows:
\begin{align}
\label{eq:mifgd} \tag{5}
U_{i+1} &= Z_{i} - \eta \mathcal{A}^\dagger \left(\mathcal{A}(Z_i Z_i^\dagger) - y\right) \cdot Z_i, \quad 
Z_{i+1} = U_{i+1} + \mu \left(U_{i+1} - U_i\right). 
\end{align}

Here, $$Z_i \in \mathbb{C}^{d\times r}$$ is a rectangular matrix (with the same dimension as $$U_i$$) which accumulates the "momentum" of the iterates $$U_i$$. $$\mu$$ is the momentum parameter that balances the weight between the previous estimate $$U_i$$ and the current estimate $$U_{i+1}.$$

While the $$\texttt{MiFGD}$$ algorithm in Eq. \eqref{eq:mifgd} looks quite similar to $$\texttt{FGD}$$ in Eq. \eqref{eq:fgd}, it complicates the convergence theory significantly. This is because the two-step momentum procedure has to be considered, on top of the fact that the objective function in Eq. \eqref{eq:factored-obj} is non-convex. We will not get into the details of the convergence thoery here; interested readers are referred to our paper.[^kim2021fast] We finish this section with an informal theorem that illustrates the convergence behavior of $$\texttt{MiFGD}$$:

**Theorem 1** ($$\texttt{MiFGD}$$ convergence rate (informal)). Assume that $$\mathcal{A}(\cdot)$$ satisfies the RIP for some constant $$0 < \delta_{2r} <1$$. Let $$y = \mathcal{A}(\rho^\star)$$ denote the set of measurements, by measuring the quantum density matrix $$\rho^\star$$. Given a good initialization point $$U_0$$, and setting step size $$\eta$$ and momentum $$\mu$$ appropriately, $$\texttt{MiFGD}$$ converges with a linear rate to a region—with radius that depends on $$O(\mu)$$—around the global solution $$\rho^\star$$. 

## Results
In this section, we review some of the experimental results. First, we obtain real quantum data from IBM's Quantum Processing Unit (QPU) by realizing two types of quantum states: $$\texttt{GHZminus}(n)$$ and $$\texttt{Hadamard}(n)$$, for $$n = 6, 8$$, where $$n$$ is the number of qubits. In quantum computing, obtaining measurements itself is not a trivial process, which we will not get into the detail in this post. Yet, we highlight that, in the following plots, we only use $$20$$% of the measurements that are information-theoretically compelete, i.e. we sample $$m = 0.2 \cdot d^2$$ measurements (recall that we are working on compressed sensing QST setting). We compare the effect of different momentum parameters in the figure below, where the accuracy of the estimated density matrix $$\widehat{\rho}$$ is measured with the true density matrix $$\rho^\star$$ in terms of the squared Frobenius norm, i.e. $$||\widehat{\rho} - \rho^\star||_F^2$$: 


![MiFGD performance on real quantum data from IBM QPU. Top-left: GHZminus(6), Top-right: GHZminus(8), Bottom-left: Hadamard(6), Bottom-right: Hadamard(8).](/assets/img/ibm-data.png)

*MiFGD performance on real quantum data from IBM QPU. Top-left: GHZminus(6), Top-right: GHZminus(8), Bottom-left: Hadamard(6), Bottom-right: Hadamard(8).*


Above figure summarizes the performance of $$\texttt{MiFGD}$$. In the legends, $$\mu^\star$$ is the momentum parameter proposed by our theory; however, it should be noted that $$\texttt{MiFGD}$$ converges with larger momentum values than $$\mu^\star$$, in particular featuring a steep dive to convergence for the largest value of $$\mu$$ we tested. Moreover, the above figure also highlights the universality of our approach: its performance is oblivious to the quantum state reconstructed, as long as it satisfies purity or it is close to a pure state. Our method does not require any additional structure assumptions in the quantum state. 

It should be noted that quantum data are inherently noisy. To highlight the level of noise existing in real quantum data, in the figure below, we also plot the performance of $$\texttt{MiFGD}$$ in the same setting using simulated quantum data from IBM's [QASM](https://github.com/Qiskit/openqasm) simulator:  

![MiFGD performance on synthetic data using IBM's QASM simulator. Top-left: GHZminus(6), Top-right: GHZminus(8), Bottom-left: Hadamard(6), Bottom-right: Hadamard(8).](/assets/img/simulator-data.png)

*MiFGD performance on synthetic data using IBM's QASM simulator. Top-left: GHZminus(6), Top-right: GHZminus(8), Bottom-left: Hadamard(6), Bottom-right: Hadamard(8).*

We see a similar trend with the result using real quantum data from IBM's QPU. However, we see that the overall accuracy of the reconstucted and the target states, $$\|\hat{\rho} - \rho^\star\|_F^2$$, is generally lower for the real quantum data--they do not reach the accuracy level of $$10^{-1}$$, which is acchieved for all cases using QASM simulator. This difference is summarized in the figure below:

![Final fidelity of MiFGD comparison using real quantum data from IBM's QPU and simulated quantum data using QASM.](/assets/img/qpu-vs-qasm.png)

*Final fidelity of MiFGD comparison using real quantum data from IBM's QPU and simulated quantum data using QASM.*

#### Performance comparison with QST methods in $\texttt{Qiskit}$
Now, we compare the performance of $$\texttt{MiFGD}$$ with [QST methods](https://qiskit.org/documentation/stubs/qiskit.ignis.verification.TomographyFitter.html) from [$$\texttt{Qiskit}$$](https://qiskit.org/), again using IBM's QASM simulator. $$\texttt{Qiskit}$$ provides two QST methods: $$i)$$ the $$\texttt{CVXPY}$$ method which relies on convex optimiztion, and $$ii)$$ the $$\texttt{lstsq}$$ which ruses least-squares fitting. Both methods solve the following full tomography problem (not compressed sensing QST problem):

\begin{align}
 \min_{\rho \in \mathbb{C}^{d \times d}}
 \quad & f(\rho) := \tfrac{1}{2} ||\mathcal{A}(\rho) - y||_2^2 \quad 
 \text{subject to}
 \quad & \rho \succeq 0, ~\texttt{Tr}(\rho) = 1.
\end{align}

We note that $$\texttt{MiFGD}$$ is not restricted to "tall" $$U$$ scenarios to encode PSD and rank constraints: even without rank constraints, one could still exploit the matrix decomposition $$\rho = UU^\dagger$$ to avoid the PSD projection, $$\rho \succeq 0$$, where $$U \in \mathbb{C}^{d \times d}$$.

We consider the following cases: $$\texttt{GHZ}(n), \texttt{Hadamard}(n),$$ and $$\texttt{Random}(n)$$ for $$n = 3, \dots, 8$$. 
The results are shown in the figure below:

![Performance comparison with Qiskit methods. All experiments are performed on a 13” Macbook Pro with 2.3 GHz Quad-Core Intel Core i7 CPU and 32 GB RAM.](/assets/img/qiskit-comparison-plot.png)

*Performance comparison with Qiskit methods. All experiments are performed on a 13” Macbook Pro with 2.3 GHz Quad-Core Intel Core i7 CPU and 32 GB RAM.*

Some notable remarks: $$i)$$ For small-scale scenarios ($$n=3, 4$$), $$\texttt{CVXPY}$$ and $$\texttt{lstsq}$$ attain almost perfect fidelity, while being comparable or faster than $$\texttt{MiFGD}$$. $$ii)$$The difference in performance becomes apparent from $$n = 6$$ and on: while $$\texttt{MiFGD}$$ attains $$98$$% fidelity in less than $$5$$ seconds, $$\texttt{CVXPY}$$ and $$\texttt{lstsq}$$ require up to hundreds of seconds to find a good solution. Finally, while $$\texttt{MiFGD}$$ gets to high-fidelity solutions in seconds for $$n = 7, 8$$, $$\texttt{CVXPY}$$ and $$\texttt{lstsq}$$ methods require more than 12 hours execution time; however, their execution never ended, since the memory usage exceeded the system's available memory.


#### Performance comparison with neural-network QST using $$\texttt{Qucumber}$$
Next, we compare the performance of $$\texttt{MiFGD}$$ compare with neural-network based QST methods, proivded by [$$\texttt{Qucumber}$$](https://qucumber.readthedocs.io/en/stable/). We consider the same quantum states as with $$\texttt{Qiskit}$$ experiments, but here we consider the case where only $$50$$% of the measurements are available. 

We report the fidelity of the reconstruction as a function of elapsed training time for $$n = 3, 4$$ in the figure below for all methods provided by $$\texttt{Qucumber}$$: $$\texttt{PRWF}, \texttt{CWF}$$, and $$\texttt{DM}$$. Note that Time (secs) on $$x$$-axis is plotted with log-scale.

![Performance comparison with Qucumber methods. All experiments are performed on a NVidia GeForce GTX 1080 TI with 11GB RAM.](/assets/img/nn-comparison-plot.png)

*Performance comparison with Qucumber methods. All experiments are performed on a NVidia GeForce GTX 1080 TI with 11GB RAM.*

We observe that for all cases, $$\texttt{Qucumber}$$ methods are orders of magnitude slower than $$\texttt{MiFGD}$$.
For the $$\texttt{Hadamard}(n)$$ and $$\texttt{Random}(n)$$, reaching reasonable fidelities is significantly slower for both $$\texttt{CWF}$$ and $$\texttt{DM},$$ while $$\texttt{PRWF}$$ hardly improves its performance throughout the training. For the $$\texttt{GHZ}$$ case, $$\texttt{CWF}$$ and $$\texttt{DM}$$ also shows *non-monotonic* behaviors: even after a few thousands of seconds, fidelities have not "stabilized", while $$\texttt{PRWF}$$ stabilizes in very low fidelities. In comparison, $$\texttt{MiFGD}$$ is several orders of magnitude faster than both $$\texttt{CWF}$$ and $$\texttt{DM}$$ and fidelity smoothly increases to comparable or higher values. What is notable is the scalability of $$\texttt{MiFGD}$$ compared to neural network approaches for higher values of $$n$$. 

To see this more clearly, in the table below, we report the final fidelities (within the $$3$$ hour time window), and reported times. We see that for many cases, $$\texttt{CWF}$$ and $$\texttt{DM}$$ methods did not complete a single iterations within $$3$$ hours. 

![Comparison with Qucumber methods.](/assets/img/nn-comparison.png)

*Comparison with Qucumber methods.*

#### The effect of parallelization in $$\texttt{MiFGD}$$
We also study the effect of paralleization in running $$\texttt{MiFGD}$$. We parallelize the iteration step across a number of processes, that can be either distributed and network connected, or sharing memory in a multicore environment. Our approach is based on [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) specification. We assign to each process a subset of the measurement labels consumed by the parallel computation. At each step, a process first computes the local gradient using a subset of measurements. These local gradients are then communicated so that they can be added up to form the full gradient, and the full gradient is shared with each worker.

In our first round of experiments, we investigate the scalability of the parallelization approach. We vary the number $$p$$ of parallel processes $$(p=1, 2, 4, 8, 16, 32, 48, 64, 80, 96)$$, collect timings for reconstructing $$\texttt{GHZ}(4)$$, $$\texttt{Random}(6)$$ and $$\texttt{GHZminus}(8)$$ states, and report speedups $$T_p/T_1$$ we gain from $$\texttt{MiFGD}$$ in the figure bloew (left panel). We observe that the benefits of parallelization are pronounced for bigger problems (here: $$n=8$$ qubits) and maximum scalability results when we use all physical cores ($48$$ in our platform).

![Effect of parallelization of MiFGD. Left: scalability of parallelization of MiFGD for different number of processors. Middle: fidelity versus time consued for different number of processors on Hadamard(10) state. Right: The effect of momentum on Hadamard(10) state with 48 processors](/assets/img/parallel.png)

*Effect of parallelization of MiFGD. Left: scalability of parallelization of MiFGD for different number of processors. Middle: fidelity versus time consued for different number of processors on Hadamard(10) state. Right: The effect of momentum on Hadamard(10) state with 48 processors.*

Further, we move to larger problems ($$n=10$$ qubits: reporting on reconstructing $$\texttt{Hadamard}(10)$$ state) and focus on the effect parallelization to achieving a given level of fidelity in reconstruction. In the middle panel of the figure above, we illustrate the fidelity as a function of the time spent in the iteration loop of $$\texttt{MiFGD}$$ for ($$p=8, 16, 32, 48, 64$$): we observe the smooth path to convergence in all $$p$$ counts which again minimizes compute time for $$p=48$$. Note that in this case we use $$10$$% of the complete measurements, and the momenutum parameter $$\mu=\frac{1}{4}$$.

Finally, in the right panel of the figure above, we fix the number of processes to $$p=48$$, in order to minimize compute time and increase the percentage of used measurements to $$20$$% of the complete measurements available for $$\texttt{Hadamard}(10)$$. We vary the momentum parameter from $$\mu=0$$ (no acceleration) to $$\mu=\frac{1}{4}$$, and confirm that we indeed get faster convergence times in the latter case while the fidelity value remains the same (i.e. coinciding upper plateau value in the plots). We can also compare with the previous fidelity versus time plot, where the same $$\mu$$ but half the measurements are consumed: more measurements translate to faster convergence times (plateau is reached roughly $$25$$% faster; compare the green line with the yellow line in the previous plot).

## Conclusion
We have introduced the $$\texttt{MiFGD}$$ algorithm for the factorized form of the low-rank QST problems. We proved that, under certain assumptions on the problem parameters, $$\texttt{MiFGD}$$ converges linearly to a neighborhood of the optimal solution, whose size depends on the momentum parameter $$\mu$$, while using acceleration motions in a non-convex setting. We demonstrate empirically, using both simulated and real data, that $$\texttt{MiFGD}$$ outperforms non-accelerated methods on both the original problem domain and the factorized space, contributing to recent efforts on testing QST algorithms in real quantum data.
These results expand on existing work in the literature illustrating the promise of factorized methods for certain low-rank matrix problems. 
Finally, we provide a publicly available implementation of our approach, compatible to the open-source software $$\texttt{Qiskit}$$, where we further exploit parallel computations in $$\texttt{MiFGD}$$ by extending its implementation to enable efficient, parallel execution over shared and distributed memory systems.




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
