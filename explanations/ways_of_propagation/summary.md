# Density Matrix Propagation Equation

The density matrix can be reconstructed from the focal series images, utilizing the way of propagation of the density matrix. Several aspects of propagation, including partial spatial coherence equation formulism and summation of pure states formulism, enables a comprehensive investigation of the density matrix propagation.

## Foundations

With the concept of Fresnel propagation analogous to the wave function but accounting for the partial coherence, the propagation mechanisms of the partial coherent electrons are interpreted in two ways: summation of pure states formulism and spatial coherence equation formulism.

#### Derivations from the summation of pure states formulism

In general, the density operator $\hat{\rho}$ can be expressed as the sum of the outer product of the pure states $|\psi_i\rangle$. 

$$
\hat{\rho}=\sum_i c_i|\psi_i\rangle \langle\psi_i|
$$

while the density matrix $\rho(k,k')=\langle k|\hat{\rho}|k'\rangle$, it can be derived that when the density matrix propagates with a distance of $\Delta z$,

$$
\begin{aligned}
\rho_{\Delta z}(k,k')&=\sum_i c_i\langle k|\psi_{i,\Delta z}\rangle\langle\psi_{i,\Delta z}|k'\rangle\\
&=\sum_{i}c_i e^{-i\pi\lambda \Delta zk^2}\langle k|\psi_i\rangle e^{i\pi\lambda \Delta zk'^2}\langle\psi_i|k'\rangle\\
&=\sum_ie^{-i\pi\lambda\Delta z(k^2-k'^2)}c_i\langle k|\psi_{i,\Delta z}\rangle\langle\psi_{i,\Delta z}|k'\rangle\\
&=e^{-i\pi\lambda\Delta z(k^2-k'^2)}\rho(k,k')
\end{aligned}
$$

#### Derivations from the spatial coherence equation formulism

In general, not only the coherence but also other factors such as the spherical aberration and the decoherence-induced envolope function. 