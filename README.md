# Regularization scheme dependence of the behavior of wavefunction

Regularization scheme dependence of the behavior of wavefunction in the presence of a delta potential is
explored in this project by regularizing the delta function to the form δ = v_0*exp(−λ*ρ**2). The project analyzes
a non-relativistic system involving two identical particles interacting with a delta potential. The energy eigen
function of this system, written in spherical coordinates , consists of a radial part R(r) that satisfies a differential
equation subject to a delta potential. Numerical simulations for this system were performed using the Numerov
method, the main objectives included determining binding energies, verifying the constancy of the regularization
parameter BE = 2.2 MeV, and analyzing the convergence point (ρc), where the numerical solution converges
with the exp(−γ ∗ ϵbρ), which is the Analytical solution for ρ → ∞ (universal behavior). This is done initially
for regularization parameters v0 = 974.443Mev, λ = 8fm**−2 and the analysis is repeated on various pairs of {v0, λ}.
The convergence point ρ_c is identified for each regularization. Results highlight the sensitivity of the
wavefunction’s behavior to regularization schemes, depicted through graphs, offering insight insights into the
quantum mechanical systems with delta potential


Numerov’s method is a standard numerical technique that leverages on numerical integration of wavefunction and
the condition that wavefunction is both continuous and differentiable at any point ρ_m.

Regularization of delta functions is a technique used to handle the issue of infinities that arise in the calculation of
physical quantities due to the mathematical nature of delta functions. One way of regularization involves modifying
the delta function in a way that removes these infinities while preserving the physical properties of the system. This
is the approach that is adopted in this project.

Delta potentials and delta functions play an important role in Quantum physics and related areas serving as
initial approximations for modeling systems such as quantum wells and nuclear interactions. This project considers
a non-relativistic, 3-d system of 2-identical particles interacting with delta potential and explores regularisations of
the form δ = v_0*exp(−λ*ρ**2).
