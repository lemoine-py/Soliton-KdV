# Project Spectral method: Solitons and the Korteweg-de Vries equation

## Introduction

The Korteweg–de Vries equation is a **second-order nonlinear partial differential equation** that models wave propagation in shallow water. It is particularly notable for its analytical solutions, which describes waves known as "solitons", stable wave packets that maintain their shape while traveling over long distances.

This project focuses on the interaction between two solitons, a phenomenon where their nonlinear nature becomes clear. Specifically, we aim to analyze and compare the numerical solution of two interacting solitons with the superposition of their individual analytical solutions. This comparison provides insights into how nonlinear interactions influence soliton dynamics and helps evaluate the accuracy of numerical methods in capturing such effects.

## The Korteweg-de Vries equation

Its "canonical" form is

$$ \frac{\partial u}{\partial t}+6 u \frac{\partial u}{\partial x}+\frac{\partial^{3} u}{\partial x^{3}}=0 $$

This equation admits an analytical solution for a single soliton wave

$$ u(x, t)=\frac{c}{2} \cosh ^{-2}\left[\frac{\sqrt{c}}{2}(x-c t-a)\right] $$

travelling at constant speed $c$, starting from $x=a$ at $t=0$.

Note that the two-solitons case does not have an analytical solution.

## NUMERICAL METHOD

The abovementioned (first) equation on a periodic domain of size L can be simulated with a split-operator or *split-step Fourier pseudo-spectral method*. Like in the case of the nonlinear Schroedinger equation, the evolution equation is given by the application of two evolution operators, one linear and one nonlinear

$$ \frac{\partial u}{\partial t}=\mathscr{L} u+\mathscr{N} u, \quad \text { with } \mathscr{L} u=-\frac{\partial^{3} u}{\partial x^{3}} \text { and } \mathscr{N} u=-6 u \frac{\partial u}{\partial x}=-3 \frac{\partial u^{2}}{\partial x} $$

One can think to create an approximated time stepping scheme by applying first one evolution operator, and then on the output of that apply the other evolution operator.

Contrary to the case of the nonlinear Schroedinger equation, here only the linear dynamics has analytical solution: in spectral space we have

$$ \frac{\partial u}{\partial t}=\mathscr{L} u=-\frac{\partial^{3} u}{\partial x^{3}} \quad \rightarrow \quad \hat{u}\_{k}(t+\Delta t)=e^{i\left(\frac{2 \pi}{L} k\right)^{3} \Delta t} \hat{u}_{k}(t) $$

However we can still use the split-step idea. In this case we first perform an update of the linear part in spectral space, and then advance the nonlinear part with a time-stepping scheme of choice. Indicating with $\mathscr{F}$ and $\mathscr{F}^{-1}$ the discrete Fourier transform and its inverse, the update operation from time $t$ to time $t+\Delta t$ is the following:

1) given $u(x, t)$, compute the Fourier transform

    $$ \hat{u}\_{k}(t) = \mathscr{F} [u(x, t)] $$

2) advance the linear part by $\Delta t$ computing the partial update in spectral space

    $$ \hat{g}\_{}(t ; \Delta t) = e^{i\left(\frac{2 \pi}{L} k\right)^{3} \Delta t} \hat{u}_{k}(t) $$

3) apply the inverse Fourier transform obtaining the partial update in physical space, compute its square, and compute the spatial derivative of that with the spectral method

$$ \begin{aligned}
g(x, t ; \Delta t) & =\mathscr{F}^{-1}\left[\hat{g}_{k}(t ; \Delta t)\right] \\
\frac{\partial g^{2}(x, t ; \Delta t)}{\partial x} & =\mathscr{F}^{-1}\left[i \frac{2 \pi}{L} k \mathscr{F}\left[g^{2}(x, t ; \Delta t)\right]\right]
\end{aligned} $$

1) apply a time stepping scheme of choice to advance the nonlinear part in physical space and obtain the fully updated solution. For example, with 1st order Euler forward

$$ u(x, t+\Delta t)=g(x, t ; \Delta t)-3 \frac{\partial g^{2}(x, t ; \Delta t)}{\partial x} \Delta t $$

Then go back to point 1 to repeat the cycle for the next timestep.

## PROJECT DESCRIPTION

1) write a code for the Korteweg-de Vries equation using the split-step method
    in periodic boundary conditions, <br> and run it to reproduce Figure 1.
    Take initial conditions

$$ u(x, 0)=\frac{c_{1}}{2} \cosh ^{-2}\left[\frac{\sqrt{c_{1}}}{2}\left(x-a_{1} L\right)\right]+\frac{c_{2}}{2} \cosh ^{-2}\left[\frac{\sqrt{c_{2}}}{2}\left(x-a_{2} L\right)\right] $$

on $x \in[0, L]$, with $L=50,\left(c_{1}=0.75, a_{1}=0.33\right)$ and $\left(c_{2}=0.4, a_{2}=0.65\right)$. In terms of numerical setup, you can discretize the domain in $\mathrm{N}=256$ grid points and
take a timestep $\Delta t=0.0004$. This initial condition corresponds to the linear combination of the analytical solutions for two separate soliton waves at $t=0$. Plot a few snapshots (or maybe make a video if you have time) to discuss what happens when the fast soliton collides with the slow soliton.

2) since the equation is nonlinear, this initial condition does not evolve like the linear combination of the analytical solutions for the two solitons for $t>0$. To visualise this, plot a space-time plot for the the linear combination of the analytical solutions for the two solitons. Discuss the comparison between the two figures. Is there a range when you can approximate the numerical solution with the linear combination of the analytical solutions for the two waves? Can you make a quantitative comparison between numerical and analytical solution limited to this range? Then take an initial condition corresponding to only one of the two solitons, and make a full quantitative comparison between numerical and analytical solution (that now is well defined).
