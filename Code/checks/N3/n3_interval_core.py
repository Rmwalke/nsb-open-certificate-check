"""Exact fixed-dyadic enclosures for the authorized finite N3 Galerkin ODE.

An interval is an integer pair (lo, hi), representing [lo/Q, hi/Q].
No binary floating-point operations participate in a certificate.  The JSON
convention is x'_k = -shell[k]*x_k - sum coefficient*x_i*x_j.

The endpoint uses variation of constants, so the dissipative diagonal is
multiplied by exp(-lambda*h), rather than widened by interval subtraction.
All Taylor remainders are universal derivative enclosures on a validated
Picard box.  This module does not run an integration on import.
"""

from fractions import Fraction
from math import isqrt


def _ceildiv(a, b):
    assert b > 0
    return -((-a) // b)


class Context:
    def __init__(self, P=224):
        self.P = int(P)
        self.Q = 1 << self.P
        self.zero = (0, 0)
        self.one = (self.Q, self.Q)
        self._exp_cache = {}

    def F(self, endpoint):
        """Convert one encoded integer endpoint to its exact rational value."""
        return Fraction(endpoint, self.Q)

    def bounds(self, x):
        return self.F(x[0]), self.F(x[1])

    def iv(self, lo, hi=None):
        lo = Fraction(lo)
        hi = lo if hi is None else Fraction(hi)
        if lo > hi:
            raise ValueError('reversed interval')
        return ((lo.numerator * self.Q) // lo.denominator,
                _ceildiv(hi.numerator * self.Q, hi.denominator))

    @staticmethod
    def add(x, y):
        return x[0] + y[0], x[1] + y[1]

    @staticmethod
    def sub(x, y):
        return x[0] - y[1], x[1] - y[0]

    @staticmethod
    def neg(x):
        return -x[1], -x[0]

    def mul(self, x, y):
        a, b = x
        c, d = y
        products = (a*c, a*d, b*c, b*d)
        return min(products) // self.Q, _ceildiv(max(products), self.Q)

    def mulrat(self, x, r):
        r = Fraction(r)
        a, b = x[0] * r.numerator, x[1] * r.numerator
        if a > b:
            a, b = b, a
        return a // r.denominator, _ceildiv(b, r.denominator)

    def divrat(self, x, r):
        r = Fraction(r)
        if not r:
            raise ZeroDivisionError
        return self.mulrat(x, 1 / r)

    def div(self, x, y):
        if y[0] <= 0 <= y[1]:
            raise ZeroDivisionError('interval denominator contains zero')
        inv = self.iv(Fraction(self.Q, y[1]), Fraction(self.Q, y[0]))
        return self.mul(x, inv)

    def sqr(self, x):
        a, b = x
        lo = 0 if a <= 0 <= b else min(a*a, b*b)
        hi = max(a*a, b*b)
        return lo // self.Q, _ceildiv(hi, self.Q)

    def sqrt(self, x):
        if x[0] < 0:
            raise ValueError('negative lower endpoint in sqrt')
        a, b = x[0]*self.Q, x[1]*self.Q
        lo, hi = isqrt(a), isqrt(b)
        if hi*hi != b:
            hi += 1
        return lo, hi

    @staticmethod
    def hull(x, y):
        return min(x[0], y[0]), max(x[1], y[1])

    @staticmethod
    def intersect(x, y):
        z = max(x[0], y[0]), min(x[1], y[1])
        if z[0] > z[1]:
            raise ArithmeticError('empty intersection of proven enclosures')
        return z

    @staticmethod
    def contains(outer, inner, strict=False):
        if strict:
            return outer[0] < inner[0] and inner[1] < outer[1]
        return outer[0] <= inner[0] and inner[1] <= outer[1]

    def _exp_point(self, x):
        """Rigorous exp of an exact Fraction, using alternating series."""
        x = Fraction(x)
        cached = self._exp_cache.get(x)
        if cached is not None:
            return cached
        if not x:
            return self.one
        a, reductions = abs(x), 0
        while a > Fraction(1, 2):
            a /= 2
            reductions += 1
        # e^(-a), a in [0,1/2]. Consecutive alternating partial sums bracket it.
        term = total = Fraction(1)
        n = 0
        tolerance = Fraction(1, 1 << (self.P + 32))
        while True:
            n += 1
            term *= -a / n
            other = total + term
            if abs(term) <= tolerance:
                lo, hi = min(total, other), max(total, other)
                break
            total = other
        if x > 0:
            lo, hi = 1/hi, 1/lo
        result = self.iv(lo, hi)
        for _ in range(reductions):
            result = self.sqr(result)
        self._exp_cache[x] = result
        return result

    def exp(self, x):
        """Monotone interval exp; valid beyond the intended range [-50,50]."""
        lo = self._exp_point(self.F(x[0]))[0]
        hi = self._exp_point(self.F(x[1]))[1]
        return lo, hi

    def exp_rat(self, x):
        return self._exp_point(Fraction(x))


class Solver:
    def __init__(self, ctx, data, h=Fraction(1, 128), order=16):
        self.ctx = ctx
        self.h = Fraction(h)
        self.order = int(order)
        self.shells = tuple(Fraction(v) for v in data['shells'])
        self.dimension = len(self.shells)
        if self.h <= 0 or self.order < 1:
            raise ValueError('positive h and order required')
        if any(lam < 0 or lam*self.h >= 1 for lam in self.shells):
            raise ValueError('weights require 0 <= lambda*h < 1')
        # Each distinct unordered pair is convolved once per order.
        grouped = {}
        for k, i, j, coefficient in data['ode_coefficients']:
            pair = min(i, j), max(i, j)
            grouped.setdefault(pair, []).append((k, -Fraction(coefficient)))
        self.pair_terms = tuple((i, j, tuple(terms))
                                for (i, j), terms in grouped.items())
        self.h_iv = ctx.iv(self.h)
        self.tau_iv = ctx.iv(0, self.h)
        self.decays = {lam: ctx.exp_rat(-lam*self.h) for lam in set(self.shells)}
        self.weights = {lam: tuple(self._weight(lam, n)
                                  for n in range(self.order+1))
                        for lam in set(self.shells)}
        self.last_picard = None
        self.last_taylor = None
        self.last_remainder_taylor = None

    def initial(self):
        return (self.ctx.one,) + (self.ctx.zero,) * (self.dimension-1)

    def _weight(self, lam, n):
        """Integral exp(-lambda*(h-s))*s^n ds, bracketed by exact series.

        Terms are h^(n+1)*n!*(-lambda*h)^k/(n+k+1)!, so consecutive
        partial sums bracket this positive integral when lambda*h < 1.
        """
        term = self.h**(n+1) / (n+1)
        if not lam:
            return self.ctx.iv(term)
        total, k = term, 0
        tolerance = Fraction(1, 1 << (self.ctx.P + 32))
        while True:
            term *= -lam*self.h / (n+k+2)
            k += 1
            other = total + term
            if abs(term) <= tolerance:
                return self.ctx.iv(min(total, other), max(total, other))
            total = other

    def nonlinear_coefficient(self, coefficients, n):
        """Normalized nth Taylor coefficient of the quadratic N(x)."""
        c = self.ctx
        out = [c.zero] * self.dimension
        for i, j, terms in self.pair_terms:
            conv = c.zero
            if i == j:
                # Symmetric convolution; the middle term is an actual square.
                for r in range((n+1)//2):
                    product = c.mul(coefficients[r][i], coefficients[n-r][i])
                    conv = c.add(conv, c.mulrat(product, 2))
                if n % 2 == 0:
                    conv = c.add(conv, c.sqr(coefficients[n//2][i]))
            else:
                for r in range(n+1):
                    conv = c.add(conv, c.mul(coefficients[r][i],
                                             coefficients[n-r][j]))
            for k, factor in terms:
                out[k] = c.add(out[k], c.mulrat(conv, factor))
        return tuple(out)

    def vector_field(self, X):
        c = self.ctx
        nonlinear = self.nonlinear_coefficient([X], 0)
        return tuple(c.sub(nonlinear[k], c.mulrat(X[k], self.shells[k]))
                     for k in range(self.dimension))

    def _taylor(self, X, degree):
        c = self.ctx
        coeffs, nonlinear_coeffs = [tuple(X)], []
        for n in range(degree):
            nonlinear = self.nonlinear_coefficient(coeffs, n)
            nonlinear_coeffs.append(nonlinear)
            coeffs.append(tuple(c.divrat(c.sub(nonlinear[k],
                               c.mulrat(coeffs[n][k], self.shells[k])), n+1)
                                for k in range(self.dimension)))
        return coeffs, nonlinear_coeffs

    def taylor(self, X, degree):
        return self._taylor(X, degree)[0]

    @staticmethod
    def _inflate(x):
        padding = max(16, (x[1]-x[0])//8 + 1)
        return x[0]-padding, x[1]+padding

    def picard_box(self, X, max_iterations=80):
        c = self.ctx
        if len(X) != self.dimension:
            raise ValueError('state dimension mismatch')
        # This is a rejection guard, never a way to accept an enclosure.
        # Reject runaway inflation before a subsequent quadratic evaluation
        # can create enormous integers.  The relative allowance also permits
        # callers whose original state has a larger physical scale.
        initial_bits = max(abs(z).bit_length() for x in X for z in x)
        bit_cap = max(c.P + 20, initial_bits + 12)

        def check_size(box):
            if any(abs(z).bit_length() > bit_cap for x in box for z in x):
                raise ArithmeticError('Picard magnitude guard exceeded; reduce h')

        self.last_picard_trace = []
        check_size(X)
        field = self.vector_field(X)
        image = tuple(c.add(x, c.mul(self.tau_iv, f)) for x, f in zip(X, field))
        B = tuple(self._inflate(c.hull(x, z)) for x, z in zip(X, image))
        for iteration in range(max_iterations):
            check_size(B)
            field = self.vector_field(B)
            image = tuple(c.add(x, c.mul(self.tau_iv, f))
                          for x, f in zip(X, field))
            failed = tuple(k for k, (b, z) in enumerate(zip(B, image))
                           if not c.contains(b, z, strict=True))
            if iteration < 20:
                self.last_picard_trace.append({
                    'iteration': iteration+1,
                    'failed': failed,
                    'widths': tuple(b[1]-b[0] for b in B),
                })
            if not failed:
                self.last_picard_iterations = iteration+1
                return B
            # Preserve every component already enclosing its Picard image.
            # Unconditional repeated inflation can itself destroy an otherwise
            # converging small box while distant nonlinear modes propagate.
            failed_set = set(failed)
            B = tuple(self._inflate(c.hull(b, z)) if k in failed_set else b
                      for k, (b, z) in enumerate(zip(B, image)))
            check_size(B)
        raise ArithmeticError('strict Picard inclusion failed; reduce h')

    def step(self, X):
        """Return (endpoint, continuous slab); never substitutes a midpoint."""
        X = tuple(X)
        c, p = self.ctx, self.order
        B = self.picard_box(X)
        coefficients, nonlinear = self._taylor(X, p)
        remainder_coefficients, remainder_nonlinear = self._taylor(B, p+1)
        endpoints, slab = [], []
        for k, lam in enumerate(self.shells):
            endpoint = c.mul(self.decays[lam], X[k])
            for n in range(p):
                endpoint = c.add(endpoint,
                                 c.mul(self.weights[lam][n], nonlinear[n][k]))
            endpoint = c.add(endpoint, c.mul(self.weights[lam][p],
                                           remainder_nonlinear[p][k]))
            # Ordinary Taylor is used only for the continuous slab, not for
            # propagation of the endpoint's dissipative linear component.
            path = remainder_coefficients[p+1][k]
            for n in range(p, -1, -1):
                path = c.add(coefficients[n][k], c.mul(self.tau_iv, path))
            path = c.intersect(path, B[k])
            endpoint = c.intersect(endpoint, path)
            endpoints.append(endpoint)
            slab.append(path)
        self.last_picard = B
        self.last_taylor = coefficients
        self.last_remainder_taylor = remainder_coefficients
        return tuple(endpoints), tuple(slab)


def arithmetic_selfcheck():
    """Only basic arithmetic checks, deliberately no ODE steps."""
    c = Context(80)
    assert c.bounds(c.mul(c.iv(-2, 3), c.iv(-4, 5))) == (Fraction(-12), Fraction(15))
    assert c.sqr(c.iv(-2, 3)) == c.iv(0, 9)
    assert c.sqrt(c.iv(4)) == c.iv(2)
    third = c.divrat(c.iv(-1), 3)
    assert c.F(third[0]) <= Fraction(-1, 3) <= c.F(third[1])
    e0, em1, ep1 = c.exp_rat(0), c.exp_rat(-1), c.exp_rat(1)
    assert e0 == c.one
    assert Fraction(1, 3) < c.F(em1[0]) <= c.F(em1[1]) < Fraction(1, 2)
    assert 2 < c.F(ep1[0]) <= c.F(ep1[1]) < 3
    assert c.contains(c.mul(em1, ep1), c.one)
    assert c.contains(c.sqrt(c.iv(2)), c.iv(Fraction(1414213562, 10**9))) is False
    return True
