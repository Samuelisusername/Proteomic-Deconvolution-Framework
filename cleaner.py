

raw = r"""\[
S_{\text{rand}} = 
\begin{bmatrix}
\big| & \big| &        & \big| \\
prot_{x_1}^{\text{sig}} & prot_{x_2}^{\text{sig}} & \cdots & prot_{x_n}^{\text{sig}} \\
\big| & \big| &        & \big|
\end{bmatrix}, \quad \forall x_i = \mathcal{C}[i], \quad prot_{x_i}^{\text{sig}}, \, prot_{x_i}^{\text{sample}} \sim \mathcal{Prot}_{x_i}
\]"""

clean = ''.join(c for c in raw if ord(c) < 128)
print(clean)