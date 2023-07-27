# PhD_admission_prob_HZDR

The direct iteration method is an iterative procedure which makes it possible to cal-
culate the eigenvalues of an operator H (to be specific, a Hamiltonian). You start with
an arbitrary vector $ψ^{(0)}$ , then you build up iterations as follow:

$$\psi^{(n)} = \frac{H \psi^{(n)}}{||H\psi^{(n-1)}||}$$
where $||Hψ^{(n−1)}||$ is the norm of the vector $Hψ^{(n−1)}$ , $n$ is the iteration number. (The
norm of a vector $a$ is $||a|| = \sqrt{\left< a\right|\left. a \right>}$.

The ground-state energy is just
$$E^{(n)} ≡ E[Ψ^{(n)}] = \frac{\left< \psi^{(n)} \right| H \left|\psi^{(n)}\right>}{\left< \psi^{(n)} \right| \left.\psi^{(n)}\right>}$$
The energy of the first excited level can be calculated in essentially the same man-
ner, but taking into account an extra condition that the sought for solution should be
orthogonal to the ground state.

Now let’s assume that the Hamilton matrix of a system in a certain basis set is de-
scribed by the formula:

$$H_{n,m}(\alpha) = \sin\left( \frac{(n+m)^{0.52}}{ln(n+m+\alpha)}\right)$$

where $1\leq n, m \leq N_{max}, 1<\alpha <2$ is a parameter.

1. Fill the Hamilton matrix for $N_{max} = 100$ and a certain $\alpha$.
2. Calculate the ground-state energy of the system using the direct iteration method.To find the lowest energy, you may need to shift the hole spectrum by adding a diagonal term to the Hamiltonian, that is $H' = H − λ ∗ I$, where $λ > 0$ is a constant (it should be a rather large number, e.g. 100), and I is the identity matrix. Then the eigenvalues should be ’shifted back’ by the same value of $λ$.
3. Plot the error (a difference between the energy for a given iteration and the energy when the convergence $10^{−6}$ or better–has been reached) as a function of the iteration number.
4. Repeat calculations for the first excited level. To do this, you need to orthogonalize the solution to the ground-state vector at each iteration.
