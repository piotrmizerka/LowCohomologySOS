# LowCohomologySOS
This repository provides rigorous proofs of the existence of a **sum of squares** decomposition for the first **cohomological Laplacian** $\Delta_1$ for finitely presented groups. More precisely, it looks for a positive **spectral gap** $\lambda>0$ such that $\Delta_1-\lambda I_n$ is a sum of squares. Here the Laplacian $\Delta_1$ is the matrix over the group ring given by the prescribed presentation of the group, $G=\langle s_1,\ldots,s_n|r_1,\ldots,r_m\rangle$. The formula is
$$
\Delta_1=d_0d_0^*+d_1^*d_1,
$$
where $d_0=\left[1-s_i\right]\in\mathbb{M}_{n,1}(\mathbb{R}G)$ and $d_1$, also known as *Jacobian*, is given by $d_1=\left[\frac{\partial r_i}{\partial x_j}\right]\in\mathbb{M}_{m,n}(\mathbb{R}G)$, where $\frac{\partial r_i}{\partial x_j}\in\mathbb{R}G$ is the $(i,j)^{\text{th}}$ **Fox derivative** (for the definition of the Fox derivatives, see the original papers of Fox, [doi:10.2307/1969736](doi:10.2307/1969736) and [doi:10.2307/1969686](doi:10.2307/1969686)). The involution $*$ is given by the composition of the matrix transposition and the standard involution on the group ring $\mathbb{R}G$.  

## Group cohomology
It has been shown by Lyndon, [doi:10.2307/1969440](https://www.jstor.org/stable/1969440?origin=crossref#metadata_info_tab_contents) , that, for any $G$-module $V$, the cohomology $H^*(G,V)$ can be computed from the following complex
$$
0\longrightarrow\mathbb{R}G\xrightarrow{\left[1-s_i\right]}(\mathbb{R}G)^n\xrightarrow{\left[\frac{\partial r_i}{\partial x_j}\right]}(\mathbb{R}G)^m\longrightarrow\cdots
$$
We focus on vanishing of the first and *reducibility* of the second cohomology with unitary coefficients (the reducibility means that the differential given by $d_1$ has a closed image for every unitary representation). By the recent work of Bader and Nowak, [doi:10.1016/j.jfa.2020.108730](https://www.sciencedirect.com/science/article/pii/S0022123620302731?via%3Dihub), these two conditions are equivalent to the existence of the positive spectral gap $\lambda>0$ such that the expression $\Delta_1-\lambda I_n$ is a **sum of squares**, that is, there exist matrices $M_1,\ldots,M_l$ such that
$$
\Delta_1-\lambda I_n=M_1^*M_1+\ldots+M_l^*M_l.
$$


## The computation for $\operatorname{SL}_3(\mathbb{Z})$

This repository was developed originally for computations concerning $\operatorname{SL}_3(\mathbb{Z})$, the special linear group $3\times 3$ matrices over integers. The computations provide a proof that there exists $\lambda\geq 0.32$ such that $\Delta_1-\lambda I_6$ is a sum of squares for a specific presentation of $\operatorname{SL}_3(\mathbb{Z})$ on six generators. The article on this subject is available at [2207.02783](https://arxiv.org/abs/2207.02783), and it contains the presentation of $\operatorname{SL}_3(\mathbb{Z})$ we used and the theory standing behind the proof of the existence of the positive spectral gap $\lambda$. The script [SL_3_Z_Delta_1.jl](./scripts/SL_3_Z_Delta_1.jl) can be used to replicate the main result.
