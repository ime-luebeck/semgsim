The original equation due to Rosenfalck, eq (16) in Farina et al.:
$$V(z [mm]) = \begin{cases}
    A z^3 e^{-z} + B & \text{if } z \gt 0 \\
    0 & \text{else}
\end{cases}
$$

The same function, taking a variable in [m]:
$$V(z [m]) = \begin{cases}
    A (1000z)^3 e^{-1000z} + B & \text{if } z \gt 0 \\
    0 & \text{else}
\end{cases}
$$

Its first derivative in the negative direction
$$ \psi(z) = \frac{d}{dz} V(-z)$$
as defined by Farina et al.:
$$
\psi(z [mm]) = -A e^{z} \cdot [3 z^2 + z^3] \cdot H(-z)
$$
where $H(z)$ denotes the Heaviside step function.

The same function but taking an argument in [m]:
$$
\psi(z [m]) = \frac{d}{dz} [-A z^3 10^9 e^{1000z} + B] \quad \text{if } z \lt 0 \text{ else } 0 \\
= -3A z^2 10^9 e^{1000z} -A z^3 10^{12} e^{1000z} \quad \text{if } z \lt 0 \text{ else } 0 \\
= -A e^{1000z} \cdot [3 z^2 10^9 + z^3 10^{12}] \cdot H(-z)
$$

The Wolfram Alpha query `integral from -inf to inf of ((- A e^(1000t) * (3*t^2*10^9 + t^3*10^12) * Heaviside(-t))*e^(-j omega t))` [^1] yields
$$
\mathcal{F}\{\psi(z [m])\}(f) = \frac{6 \cdot 10^9 A j \omega}{(j\omega- 1000)^4} \\
= \frac{12j \cdot 10^9 A \pi f}{(2 j \pi f- 1000)^4}
$$
which is the quantity called $\Psi$ in the Farina paper.

Let's try to derive the same thing analytically.

We know that (for a non-unitary Fourier transform, cf. Fourier-transforms.Rmd)
$$
x^n f(x) \mapsto i^n \frac{d^n \hat{f}(\omega)}{d\omega^n}.
$$

Our expression for $\psi$ above can be decomposed as
$$
\psi(z[m]) = -A e^{1000z} \cdot [3 z^2 10^9 + z^3 10^{12}] \cdot H(-z) \\
= - 3A z^2 e^{1000z} 10^9 H(-z) - A z^3 e^{1000z} 10^{12} H(-z).
$$

...

[^1]: Notice that that query is *not* the same as simply calling `fourier transform of (- A e^(1000t) * (3*t^2*10^9 + t^3*10^12) * Heaviside(-t)`, because that will use another type of Fourier transform! (Probably a unitary one?) The difference should be a factor of $\sqrt{2 \pi}$.