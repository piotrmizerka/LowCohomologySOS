
# LowCohomologySOS

This repository provides **sum of squares** witnesses for the existence of the **spectral gap** for the first **cohomological Laplacian** $\Delta_1$ for finitely presented groups.

More precisely we are looking for a positive $\lambda$ such that $\Delta_1-\lambda I_n$ admits a sum of squares decomposition, which proves that $\lambda$ is a lower bound on the spectral gap of $\Delta_1$.

For a finitely presented group $G=\langle s_1,\ldots,s_n|r_1,\ldots,r_m\rangle$ the first Laplacian $\Delta_1$ is given by the formula

$$\Delta_1=d_0d_0^*+d_1^*d_1\in\mathbb M_{n,n}(\mathbb ZG)\subseteq M_{n,n}(\mathbb RG),$$

where $d_0=\left[1-s_i\right]\in\mathbb M_{n,1}(\mathbb ZG)$ and $d_1$, known as *Jacobian*, is given by the $m\times n$ matrix of *Fox derivatives* of the relations (for the definition see the original papers of Fox, [doi:10.2307/1969736](https://www.jstor.org/stable/1969736#metadata_info_tab_contents) and [doi:10.2307/1969686](https://www.jstor.org/stable/1969686#metadata_info_tab_contents)). The involution $*$ is given by the composition of the matrix transposition and the standard involution on the group ring $\mathbb RG$.  

## Group cohomology

It has been shown by Lyndon in [doi:10.2307/1969440](https://www.jstor.org/stable/1969440) that for any $G$-module $V$ the cohomology $H^*(G,V)$ can be computed from the following complex

$$0\longrightarrow\mathbb ZG\xrightarrow{d_0}(\mathbb ZG)^n\xrightarrow{d_1}(\mathbb ZG)^m\longrightarrow\cdots$$

We focus on the vanishing of the first and the *reducibility* of the second cohomology with unitary coefficients (here reducibility means that the the image of differential $d_1$ is closed for every unitary $G$-module $V$).
By the recent work of Bader and Nowak [doi:10.1016/j.jfa.2020.108730](https://www.sciencedirect.com/science/article/pii/S0022123620302731) these two conditions are equivalent to the existence of a positive $\lambda>0$ such that $\Delta_1-\lambda I_n$ admits a sum of squares decomposition, that is, there exist matrices $M_1,\ldots,M_l$ such that

$$
\Delta_1-\lambda I_n=M_1^*M_1+\ldots+M_l^*M_l.
$$

## The computation for $\operatorname{SL}_3(\mathbb{Z})$

## Replication for [2207.02783](https://arxiv.org/abs/2207.02783)
### Note: replication for [2207.02783](https://arxiv.org/abs/2207.02783) has been moved to a separate branch [2207.02783](https://github.com/piotrmizerka/LowCohomologySOS/tree/2207.02783).
Script [SL_3_Z_Delta_1.jl](./scripts/SL_3_Z_Delta_1.jl) provides a proof of the existence of $\lambda\geq 0.32$ such that $\Delta_1-\lambda I_6$ is a sum of squares for the Steinberg presentation of $\operatorname{SL}_3(\mathbb{Z})$ on six generators (defined in Section 2 of [2207.02783](https://arxiv.org/abs/2207.02783)).

The running time of the script will be approximately 3 hours on a standard laptop computer.

## Citing

If you use any code from this repository, or you find reading through the code enlightening please cite [2207.02783](https://arxiv.org/abs/2207.02783) as

```bash
@misc{https://doi.org/10.48550/arxiv.2207.02783,
  doi = {10.48550/ARXIV.2207.02783},  
  url = {https://arxiv.org/abs/2207.02783},  
  author = {Kaluba, Marek and Mizerka, Piotr and Nowak, Piotr W.},  
  keywords = {Group Theory (math.GR), Operator Algebras (math.OA), FOS: Mathematics, FOS: Mathematics},  
  title = {Spectral gap for the cohomological Laplacian of $\operatorname{SL}_3(\mathbb{Z})$},
  publisher = {arXiv},  
  year = {2022},  
  copyright = {arXiv.org perpetual, non-exclusive license}
}
