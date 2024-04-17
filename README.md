
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

# Replication details for [2404.10287](https://arxiv.org/abs/2404.10287)

For the computations we used julia in version `1.8.3` but in principle any later version should work.

## Obtaining code
To obtain the code first clone this repository via
```bash
git clone https://github.com/piotrmizerka/LowCohomologySOS.git
```
and checkout to the correct branch
```bash
cd LowCohomologySOS
git checkout sln_induction
```

## Setting up the environment
First, run julia in `LowCohomologySOS` folder
```bash
julia --project=.
```
Next, to set up the proper environment for the replication run in julia REPL
```julia
julia> using Pkg; Pkg.instantiate()
```
This command installs and precompiles, if needed, all the necessary dependences,
so it may take a while.
Note that this step needs to be executed only once per installation.

## Running actual replication
We wish to prove that for the Steinberg presentation of $\text{SL}_3(\mathbb{Z})$
on six generators (as defined in Section 3 of [2404.10287](https://arxiv.org/abs/2404.10287))
$\text{Adj}_3-\lambda I_6$ is a sum of squares for some $\lambda\geq 0.217$.

We provide a script which performs the necessary optimization to find such sum of squares decomposition.

As before the following command needs to be executed in `LowCohomologySOS` folder:
```bash
julia --project=. ./scripts/SL_3_Z_adj.jl
```

The running time of the script will be approximately `2` hours on a standard laptop computer.

Instead of running the whole computation, one can use the precomputed solution instead. In order to run the script providing rigorous mathematical proof (see the Section 3.2 of [2207.02783](https://arxiv.org/abs/2207.02783)) that $\text{Adj}_3-0.217 I_6$ is a sum of squares, execute the following command in `LowCohomologySOS` folder:
```bash
julia --project=. ./scripts/sl3_adj_precom/SL_3_Z_adj_cert.jl
```
The running time of the script will be approximately `2` minutes on a standard laptop computer.