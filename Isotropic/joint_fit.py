"""Joint binned Poisson fit used by ``fit_isotropic.py fit``.

Model for energy bin k and pixel p (all arrays are predicted COUNTS):

    mu[k, p] = fixed[k, p] + a_k * T_iso[k, p] + sum_j A_j * T_j[k, p]

* ``a_k``  : one isotropic amplitude per energy bin (ne parameters)
* ``A_j``  : one energy-independent amplitude per free diffuse template
* ``fixed``: everything held constant (catalog sources, fixed diffuse parts)

All amplitudes are constrained to be >= 0. Because the free diffuse amplitudes
are shared across energies, all bins are fitted simultaneously.

No instrument-response calculation takes place here: the templates must
already include exposure, PSF and (for non-isotropic components) energy
dispersion. The isotropic template must have energy dispersion disabled.
"""

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import xlogy


class JointPoisson:
    """Joint Poisson likelihood for per-bin isotropic + global diffuse amplitudes.

    Parameters
    ----------
    counts, fixed : ndarray, shape (n_energy, n_pixel)
        Observed counts and fixed model counts (retained pixels only).
    templates : ndarray, shape (1 + n_free, n_energy, n_pixel)
        ``templates[0]`` is the isotropic template; the rest are the free
        diffuse templates, each at unit amplitude.

    The parameter vector is ``[a_0 ... a_{ne-1}, A_0 ... A_{n_free-1}]``.
    """

    def __init__(self, counts, fixed, templates):
        self.n, self.fixed, self.t = [np.asarray(v, float) for v in (counts, fixed, templates)]
        if (
            self.n.ndim != 2
            or self.fixed.shape != self.n.shape
            or self.t.ndim != 3
            or self.t.shape[1:] != self.n.shape
            or len(self.t) < 1
        ):
            raise ValueError(
                "Expected counts/fixed [energy,pixel], templates [component,energy,pixel]."
            )
        if any(np.any(~np.isfinite(v)) or np.any(v < 0) for v in (self.n, self.fixed, self.t)):
            raise ValueError("Invalid counts/templates.")
        if np.any((self.fixed + self.t.sum(axis=0) == 0) & (self.n > 0)):
            raise ValueError("Counts without model support.")
        self.ne = self.n.shape[0]
        self.ng = len(self.t) - 1
        self.npar = self.ne + self.ng
        self.scale = max(float(self.n.sum()), 1.0)
        if np.any(self.t[0].sum(axis=1) <= 0) or np.any(self.t[1:].sum(axis=(1, 2)) <= 0):
            raise ValueError("Fitted template without support.")

    def model(self, a):
        """Predicted counts for parameter vector ``a``."""
        return (
            self.fixed
            + a[: self.ne, None] * self.t[0]
            + np.einsum("j,jkp->kp", a[self.ne :], self.t[1:])
        )

    def objective(self, a):
        """Return (half Poisson deviance, gradient) at ``a``.

        Half deviance = sum[mu - n + n ln(n/mu)] equals the negative log
        likelihood up to a constant, so differences are delta(-ln L). A
        model that is zero where counts are positive returns ``inf``.
        """
        mu = self.model(a)
        if np.any((mu <= 0) & (self.n > 0)):
            return np.inf, np.full(self.npar, -1e100)
        # Half the Poisson deviance differs from NLL by a constant.
        value = np.sum(mu - self.n + xlogy(self.n, self.n) - xlogy(self.n, mu))
        r = 1.0 - np.divide(self.n, mu, out=np.zeros_like(mu), where=mu > 0)
        g = np.r_[(r * self.t[0]).sum(axis=1), np.einsum("kp,jkp->j", r, self.t[1:])]
        return float(value), g

    def covariance(self, a):
        """Inverse of the Hessian sum_p n/mu^2 * t_i t_j at ``a``.

        This is the usual Gaussian (curvature) covariance of the amplitudes.
        Isotropic bins are coupled only through the shared diffuse amplitudes.
        Returns an all-NaN matrix if the Hessian is ill conditioned (cond > 1e12),
        and is not meaningful for parameters fitted at the boundary (0).
        """
        mu = self.model(a)
        w = np.divide(self.n, mu**2, out=np.zeros_like(mu), where=mu > 0)
        h = np.zeros((self.npar, self.npar))
        h[np.arange(self.ne), np.arange(self.ne)] = (w * self.t[0] ** 2).sum(axis=1)
        cross = np.einsum("kp,kp,jkp->kj", w, self.t[0], self.t[1:])
        h[: self.ne, self.ne :] = cross
        h[self.ne :, : self.ne] = cross.T
        h[self.ne :, self.ne :] = np.einsum("kp,ikp,jkp->ij", w, self.t[1:], self.t[1:])
        if np.linalg.cond(h) > 1e12:
            return np.full(h.shape, np.nan)
        return np.linalg.inv(h)

    def solve(self, start=None, fixed_parameter=None):
        """Minimize the objective with L-BFGS-B and nonnegative bounds.

        Two starting points are tried (``start`` and all ones); a solution is
        accepted only if the optimizer converged and the projected gradient
        passes a KKT check. The best accepted solution is returned.

        Parameters
        ----------
        start : array-like, optional
            Initial parameters (default: all ones).
        fixed_parameter : (index, value), optional
            Hold one parameter at ``value`` (used for profiling).

        Returns
        -------
        a : ndarray
            Best-fit parameters.
        value : float
            Half deviance at the minimum.
        """
        start = np.ones(self.npar) if start is None else np.asarray(start, float).copy()
        free = np.ones(self.npar, bool)
        if fixed_parameter is not None:
            i, v = fixed_parameter
            free[i] = False
            start[i] = v

        def expand(x):
            a = start.copy()
            a[free] = x
            return a

        def fg(x):
            # Scaled by total counts so the stopping tolerances are dimensionless.
            v, g = self.objective(expand(x))
            return v / self.scale, g[free] / self.scale

        if not np.any(free):
            return start, self.objective(start)[0]
        solutions = []
        for x in [start[free], np.ones(free.sum())]:
            r = minimize(
                fg,
                x,
                jac=True,
                bounds=[(0, None)] * len(x),
                method="L-BFGS-B",
                options={
                    "ftol": 1e-15,
                    "gtol": 2e-11,
                    "maxiter": 5000,
                    "maxls": 60,
                    "maxcor": min(30, len(x)),
                },
            )
            g = fg(r.x)[1]
            # Projected gradient: ignore components pushing against the a >= 0 bound.
            pg = np.where((r.x < 1e-8) & (g > 0), 0.0, g)
            if r.success and np.isfinite(r.fun) and np.max(np.abs(pg)) < 2e-7:
                solutions.append(r)
        if not solutions:
            raise RuntimeError("Joint optimizer failed convergence/KKT checks.")
        r = min(solutions, key=lambda r: r.fun)
        a = expand(r.x)
        return a, self.objective(a)[0]

    def profile(self, index, best, minimum):
        """Profile-likelihood summary for parameter ``index``.

        All other parameters are re-fitted at each trial value.

        Returns
        -------
        lo, hi : float
            68% interval (delta -lnL = 0.5). ``lo`` is 0 if the interval
            reaches zero.
        ul95 : float
            One-sided 95% upper limit (delta -lnL = 2.71/2 = 1.3528).
        ts : float
            Test statistic, 2 * [-lnL(parameter = 0) - -lnL(best)].
        """

        def delta(v):
            # Profiled -lnL(v) - -lnL(best). First check that setting the
            # parameter to v cannot make mu = 0 where counts are positive.
            # A zero component can make a positive-count pixel impossible.
            a = best.copy()
            a[index] = v
            support = self.fixed.copy()
            if index < self.ne:
                support += self.t[0]
                support[index] = self.fixed[index] + v * self.t[0, index]
                support += self.t[1:].sum(axis=0)
            else:
                support += self.t[0] + self.t[1:].sum(axis=0)
                support += (v - 1.0) * self.t[index - self.ne + 1]
            if np.any((support <= 0) & (self.n > 0)):
                return np.inf
            _, f = self.solve(best, (index, v))
            return f - minimum

        a = best[index]
        d0 = delta(0.0)  # used for TS and to decide whether lo = 0
        if d0 < -1e-3:
            raise RuntimeError("Profile found a better minimum; global fit must be refined.")

        def upper(level):
            # Expand the bracket until delta exceeds `level`, then root-find.
            hi = max(1.0, 2 * a + 0.1)
            for _ in range(60):
                if delta(hi) >= level:
                    return brentq(lambda v: delta(v) - level, a, hi, xtol=1e-7)
                hi *= 2
            raise RuntimeError("Unbounded profile interval / non-identifiable component.")

        lo = 0.0 if d0 <= 0.5 or a < 1e-10 else brentq(lambda v: delta(v) - 0.5, 0.0, a, xtol=1e-7)
        return lo, upper(0.5), upper(1.352771727), max(0.0, 2 * d0)
