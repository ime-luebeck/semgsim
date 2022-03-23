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

The Wolfram Alpha query `fourier transform of - A e^(t) * (3*t^2 + t^3) * Heaviside(-t)` yields
$$
\mathcal{F}\{\psi(z [mm])\} = \frac{-3i \sqrt{2/\pi} A \omega}{(\omega - i)^4}
$$

Knowing that
$$
\mathcal{F}\{f(\alpha t)\}(\omega) = \frac{1}{\lvert \alpha \rvert} \mathcal{F}\{f(t)\}(\omega/\alpha),
$$
we obtain
$$
\mathcal{F}\{\psi(z [m])\}(\omega) = \frac{1}{1000}\mathcal{F}\{\psi(z [mm])\}(\omega/1000)=  \\
= \frac{-3i \cdot 10^{-3} \sqrt{2/\pi} A \omega \cdot 10^{-3}}{(\omega \cdot 10^{-3}- i)^4} \\
= \frac{-3i \cdot 10^{-3} \sqrt{2/\pi} A \omega \cdot 10^{-3}}{(\pi f- 500i)^4} \cdot 500^4 \\
= \frac{-6i \cdot 10^{-6} \sqrt{2\pi} A f \cdot 500^4}{(\pi f- 500i)^4} \\
= \frac{-375000i \sqrt{2\pi} A f}{(\pi f- 500i)^4}
$$
which is what is called $\Psi$ in the Farina paper.