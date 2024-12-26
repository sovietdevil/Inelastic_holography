# Review on density matrix transform-based holography
The articles below demonstrates the way of propagation for the partial coherent electrons.
- *Falk Rödern, Axel Lubk, Transfer and reconstruction of the density matrix in off-axis electron holography, Ultramicroscopy 146 103-116, (2014)*
## Concepts
**TCC**: High-resolution images obtained from elastically scattered electrons are modulated by the combined influence of coherent lens aberrations and the  partially coherent source summarized in the transmission cross-coefficient (TCC).

**Scattering-induced decoherence**: Different inelastic
scattering processes [3–9] and statistical fluctuations within the object [10], within the energy and momentum distribution of the electron beam or within the surrounding environment [11–13] also change the coherence properties of the beam.

**MDFF**: the dynamic form factor was later generalized to the mixed dynamic form factor (MDFF) [15] incorporating also the off-diagonals (i.e. coherence).
- First Born approximation: MDFF is proportional to the cross-spectral density at the back focal plane.

**Mutual coherence function**: inverse temporal Fourier transformation of the cross-spectral density is called mutual coherence function in turn [17, 18]

## Core equations
#### Density matrix of the probe electron
Flux of the pure state
$$
\mathbf{j}=\frac{\hbar}{2 i m}\left(\psi_0^* \nabla \psi_0-\psi_0 \nabla \psi_0^*\right) \approx \frac{\hbar \mathbf{k}_0}{m}\left|\psi_0\right|^2
$$
Probability density of the electron beam
$$
\rho_s(\mathbf{r})=\sum_{i j} \psi_i(\mathbf{r}) \psi_j^*(\mathbf{r})\left\langle\tau_j \mid \tau_i\right\rangle=\sum_i\left|\psi_i(\mathbf{r})\right|^2
$$
#### Electron source
In Fourier space:
$$
\rho_b\left(\mathbf{k}_{\perp}, \mathbf{k}_{\perp}^{\prime}, E\right)=f_{\mathrm{C}}(E) f_{\mathrm{S}}\left(\mathbf{k}_{\perp}\right) \delta\left(\mathbf{k}_{\perp}-\mathbf{k}_{\perp}^{\prime}\right)
$$
In real space:
$$
\rho_b\left(\mathbf{r}, \mathbf{r}^{\prime}, E\right)=\mathrm{FT}^{-1}\left[\rho_b\left(\mathbf{k}_{\perp}, \mathbf{k}_{\perp}^{\prime}, E\right)\right]=f_{\mathrm{C}}(E) \tilde{f}_{\mathrm{S}}\left(\mathbf{r}-\mathbf{r}^{\prime}\right)
$$

## Statements
- van Hove showed in
1954 that the dynamic form factor relates the diffracted intensity
and density-density correlation in the object [14]
## References
#### Inelstic scattering processes
[3] F.J. García de Abajo, Optical excitations in electron microscopy, Rev. Mod. Phys.
82 (2010) 209–275.

[4] A. Howie, Inelastic scattering of electrons by crystals: I. The theory of small-
angle inelastic scattering, Proc. R. Soc. Lond. A 271 (1963) 268–287.

[5] L. Reimer, R. Rennekamp, Imaging and recording of multiple scattering effects
by angular resolved electron energy loss spectroscopy, Ultramicroscopy 28
(1989) 258–265.

[6] Y.Y. Wang, S.C. Cheng, V.P. Dravid, F.C. Zhang, Momentum-transfer resolved
electron energy loss spectroscopy of solids: problems, solutions and applica-
tions, Ultramicroscopy 59 (1995) 109–119.

[7] L. Gu, V. Srot, W. Sigle, C. Koch, P. van Aken, F. Scholz, S.B. Thapa, C. Kirchner,
M. Jetter, M. Rühle, Band-gap measurements of direct and indirect semicon-
ductors using monochromated electrons, Phys. Rev. B 75 (2007) 195214.

[8] H. Kohl, H. Rose, Theory of image formation by inelastically scattered electrons
in the electron microscope, Adv. Electron. Electron Phys. 65 (1985) 173–227.

[9] S.L. Dudarev, L.-M. Peng, M.J. Whelan, Correlations in space and time and
dynamical diffraction of high-energy electrons by crystals, Phys. Rev. B
48 (1993) 13408.
F. Röder, A. Lubk / Ultramicroscopy 146 (2014) 103–116 115
#### Statistical fluctuations
[10] A. Rother (Lubk), T. Gemming, H. Lichte, The statistics of the thermal motion of
the atoms during imaging process in transmission electron microscopy and
related techniques, Ultramicroscopy 109 (2009) 139–146.
#### Surrounding environment effect
[11] F.J. García de Abajo, Optical emission from the interaction of fast electrons
with metallic films containing a circular aperture: a study of radiative
decoherence of fast electrons, Phys. Rev. Lett. 102 (2009) 237401.

[12] A. Howie, Mechanisms of decoherence in electron microscopy, Ultramicro-
scopy 111 (2011) 761–767.

[13] S. Uhlemann, H. Müller, P. Hartel, J. Zach, M. Haider, Thermal magnetic field
noise limits resolution in transmission electron microscopy, Phys. Rev. Lett.
111 (2013) 046101.
