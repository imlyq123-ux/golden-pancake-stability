// Independent arbitrary-precision Gragg-Bulirsch-Stoer (GBS) extrapolation integrator
// for the planar Newtonian three-body problem with masses (1, mu, 1) and, optionally,
// its variational equations (2 shooting columns or the full 12x12 monodromy).
//
// This file deliberately shares no numerical code with taylor_mpfr.cpp: it uses only
// right-hand-side evaluations (forces and the analytic force Jacobian), Gragg's
// modified midpoint rule with smoothing, polynomial (Aitken-Neville) extrapolation in
// h^2 over the step-number sequence n_j = 2j, and an embedded extrapolation error
// estimate for step-size control. Every arithmetic operation is MPFR, rounded to nearest.
//
// Usage: gbs_mpfr digits tolerance stages mu u v T columns output stride
//   digits    decimal digits printed; working precision is digits+20 decimal digits
//   tolerance accepted when max_i |T_kk - T_k,k-1| <= tolerance * (1 + |y_i|)
//   stages    k, number of extrapolation stages (order 2k of the accepted value)
//   columns   0, 2 (d/du, d/dv of the Euler-section initial state) or 12 (monodromy)
//   output    CSV: t, 12 state components, then each variational column contiguously
#include <cstdio>   // must precede mpfr.h so that mpfr_fprintf is declared
#include <mpfr.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

static mpfr_prec_t PREC = 256;
static const mpfr_rnd_t RND = MPFR_RNDN;

struct Num {
    mpfr_t v;
    Num() { mpfr_init2(v, PREC); mpfr_set_zero(v, 1); }
    Num(const Num& o) { mpfr_init2(v, PREC); mpfr_set(v, o.v, RND); }
    Num& operator=(const Num& o) { mpfr_set(v, o.v, RND); return *this; }
    ~Num() { mpfr_clear(v); }
};
using Vec = std::vector<Num>;

struct System {
    int cols = 0, n = 12;
    Num m[3];
    // scratch
    Num dx, dy, r2, r, inv3, inv5, b00, b01, b11, t1, t2, eqx, eqy, ax, ay;

    void rhs(const Vec& y, Vec& f) {
        // positions derivative = velocities
        for (int a = 0; a < 6; a++) mpfr_set(f[a].v, y[6 + a].v, RND);
        for (int a = 6; a < 12; a++) mpfr_set_zero(f[a].v, 1);
        for (int c = 0; c < cols; c++) {
            int base = 12 + 12 * c;
            for (int a = 0; a < 6; a++) mpfr_set(f[base + a].v, y[base + 6 + a].v, RND);
            for (int a = 6; a < 12; a++) mpfr_set_zero(f[base + a].v, 1);
        }
        static const int I[3] = {0, 0, 1}, J[3] = {1, 2, 2};
        for (int p = 0; p < 3; p++) {
            int i = I[p], j = J[p];
            mpfr_sub(dx.v, y[2 * j].v, y[2 * i].v, RND);
            mpfr_sub(dy.v, y[2 * j + 1].v, y[2 * i + 1].v, RND);
            mpfr_mul(r2.v, dx.v, dx.v, RND);
            mpfr_fma(r2.v, dy.v, dy.v, r2.v, RND);
            mpfr_sqrt(r.v, r2.v, RND);
            mpfr_mul(inv3.v, r2.v, r.v, RND);
            mpfr_ui_div(inv3.v, 1, inv3.v, RND);  // |d|^-3
            // acceleration contributions
            mpfr_mul(ax.v, dx.v, inv3.v, RND);
            mpfr_mul(ay.v, dy.v, inv3.v, RND);
            mpfr_fma(f[6 + 2 * i].v, m[j].v, ax.v, f[6 + 2 * i].v, RND);
            mpfr_fma(f[7 + 2 * i].v, m[j].v, ay.v, f[7 + 2 * i].v, RND);
            mpfr_mul(t1.v, m[i].v, ax.v, RND); mpfr_sub(f[6 + 2 * j].v, f[6 + 2 * j].v, t1.v, RND);
            mpfr_mul(t1.v, m[i].v, ay.v, RND); mpfr_sub(f[7 + 2 * j].v, f[7 + 2 * j].v, t1.v, RND);
            if (!cols) continue;
            // B = |d|^-3 I - 3 |d|^-5 d d^T  (derivative of d/|d|^3 with respect to d)
            mpfr_div(inv5.v, inv3.v, r2.v, RND);
            mpfr_mul_ui(inv5.v, inv5.v, 3, RND);
            mpfr_mul(t1.v, dx.v, dx.v, RND); mpfr_mul(t1.v, t1.v, inv5.v, RND); mpfr_sub(b00.v, inv3.v, t1.v, RND);
            mpfr_mul(t1.v, dy.v, dy.v, RND); mpfr_mul(t1.v, t1.v, inv5.v, RND); mpfr_sub(b11.v, inv3.v, t1.v, RND);
            mpfr_mul(b01.v, dx.v, dy.v, RND); mpfr_mul(b01.v, b01.v, inv5.v, RND); mpfr_neg(b01.v, b01.v, RND);
            for (int c = 0; c < cols; c++) {
                int base = 12 + 12 * c;
                // e = dq_j - dq_i ; w = B e
                mpfr_sub(eqx.v, y[base + 2 * j].v, y[base + 2 * i].v, RND);
                mpfr_sub(eqy.v, y[base + 2 * j + 1].v, y[base + 2 * i + 1].v, RND);
                mpfr_mul(ax.v, b00.v, eqx.v, RND); mpfr_fma(ax.v, b01.v, eqy.v, ax.v, RND);
                mpfr_mul(ay.v, b01.v, eqx.v, RND); mpfr_fma(ay.v, b11.v, eqy.v, ay.v, RND);
                mpfr_fma(f[base + 6 + 2 * i].v, m[j].v, ax.v, f[base + 6 + 2 * i].v, RND);
                mpfr_fma(f[base + 7 + 2 * i].v, m[j].v, ay.v, f[base + 7 + 2 * i].v, RND);
                mpfr_mul(t2.v, m[i].v, ax.v, RND); mpfr_sub(f[base + 6 + 2 * j].v, f[base + 6 + 2 * j].v, t2.v, RND);
                mpfr_mul(t2.v, m[i].v, ay.v, RND); mpfr_sub(f[base + 7 + 2 * j].v, f[base + 7 + 2 * j].v, t2.v, RND);
            }
        }
    }
};

static void axpy(Vec& out, const Vec& x, const Num& h, const Vec& f) {  // out = x + h f
    for (size_t i = 0; i < out.size(); i++) mpfr_fma(out[i].v, h.v, f[i].v, x[i].v, RND);
}

int main(int argc, char** argv) {
    try {
        if (argc != 11) throw std::runtime_error("Usage: gbs_mpfr digits tolerance stages mu u v T columns output stride");
        int digits = std::atoi(argv[1]);
        int k = std::atoi(argv[3]);
        int cols = std::atoi(argv[8]);
        int stride = std::atoi(argv[10]);
        if (digits < 20 || digits > 400) throw std::runtime_error("digits must be in 20..400");
        if (k < 4 || k > 60) throw std::runtime_error("stages must be in 4..60");
        if (cols != 0 && cols != 2 && cols != 12) throw std::runtime_error("columns must be 0, 2 or 12");
        if (stride < 1) throw std::runtime_error("stride must be positive");
        PREC = (mpfr_prec_t)std::ceil((digits + 20) * std::log2(10.0));

        System sys; sys.cols = cols; sys.n = 12 * (1 + cols);
        const int N = sys.n;
        Num tol, mu, u, v, T, t, H, Hmax, two, half;
        auto parse = [](Num& x, const char* s) {
            if (mpfr_set_str(x.v, s, 10, RND) != 0) throw std::runtime_error(std::string("bad number: ") + s);
        };
        parse(tol, argv[2]); parse(mu, argv[4]); parse(u, argv[5]); parse(v, argv[6]); parse(T, argv[7]);
        if (mpfr_sgn(tol.v) <= 0 || mpfr_sgn(mu.v) <= 0 || mpfr_sgn(T.v) <= 0) throw std::runtime_error("tolerance, mu and T must be positive");
        mpfr_set_ui(sys.m[0].v, 1, RND); mpfr_set(sys.m[1].v, mu.v, RND); mpfr_set_ui(sys.m[2].v, 1, RND);
        mpfr_set_ui(two.v, 2, RND); mpfr_set_d(half.v, 0.5, RND);

        // Euler-section initial state; v2 = -2 (u, v) / mu formed at working precision.
        Vec y(N);
        mpfr_set_si(y[0].v, -1, RND); mpfr_set_si(y[4].v, 1, RND);
        mpfr_set(y[6].v, u.v, RND); mpfr_set(y[7].v, v.v, RND);
        mpfr_mul_si(y[8].v, u.v, -2, RND); mpfr_div(y[8].v, y[8].v, mu.v, RND);
        mpfr_mul_si(y[9].v, v.v, -2, RND); mpfr_div(y[9].v, y[9].v, mu.v, RND);
        mpfr_set(y[10].v, u.v, RND); mpfr_set(y[11].v, v.v, RND);
        if (cols == 12) for (int c = 0; c < 12; c++) mpfr_set_ui(y[12 + 12 * c + c].v, 1, RND);
        if (cols == 2) {
            Num m2; mpfr_si_div(m2.v, -2, mu.v, RND);
            for (int c = 0; c < 2; c++) {  // d(state)/du (c=0), d(state)/dv (c=1)
                mpfr_set_ui(y[12 + 12 * c + 6 + c].v, 1, RND);
                mpfr_set(y[12 + 12 * c + 8 + c].v, m2.v, RND);
                mpfr_set_ui(y[12 + 12 * c + 10 + c].v, 1, RND);
            }
        }

        FILE* out = std::fopen(argv[9], "w");
        if (!out) throw std::runtime_error("cannot open output");
        std::fprintf(out, "# GBS MPFR %s digits=%d bits=%ld stages=%d tol=%s columns=%d sequence=2j smoothing=gragg\n",
                     mpfr_get_version(), digits, (long)PREC, k, argv[2], cols);
        auto emit = [&](const Num& time, const Vec& state) {
            mpfr_fprintf(out, "%.*Re", digits, time.v);
            for (int i = 0; i < N; i++) mpfr_fprintf(out, ",%.*Re", digits, state[i].v);
            std::fputc('\n', out);
        };
        emit(t, y);

        // Workspace.
        Vec f0(N), fz(N), zprev(N), z(N), znext(N), scale(N);
        std::vector<Vec> row(k, Vec(N)), prev(k, Vec(N));
        std::vector<int> nseq(k);
        for (int j = 0; j < k; j++) nseq[j] = 2 * (j + 1);
        // Neville coefficients 1/((n_j/n_{j-l})^2 - 1) in MPFR.
        std::vector<std::vector<Num>> coef(k, std::vector<Num>(k));
        for (int j = 0; j < k; j++)
            for (int l = 1; l <= j; l++) {
                Num a;
                mpfr_set_ui(a.v, (unsigned long)nseq[j], RND); mpfr_div_ui(a.v, a.v, (unsigned long)nseq[j - l], RND);
                mpfr_sqr(a.v, a.v, RND); mpfr_sub_ui(a.v, a.v, 1, RND);
                mpfr_ui_div(coef[j][l].v, 1, a.v, RND);
            }

        mpfr_set_d(H.v, 1e-3, RND);
        mpfr_set_d(Hmax.v, 0.25, RND);
        long steps = 0, rejected = 0, evals = 0;
        Num h, h2, tmp, err, ratio, Hnew, remaining;
        const double expo = 1.0 / (2.0 * k - 1.0);
        while (mpfr_cmp(t.v, T.v) < 0) {
            if (steps > 2000000) throw std::runtime_error("step limit reached");
            bool last = false;
            mpfr_sub(remaining.v, T.v, t.v, RND);
            if (mpfr_cmp(H.v, remaining.v) >= 0) { mpfr_set(H.v, remaining.v, RND); last = true; }
            sys.rhs(y, f0); evals++;
            for (int i = 0; i < N; i++) { mpfr_abs(scale[i].v, y[i].v, RND); mpfr_add_ui(scale[i].v, scale[i].v, 1, RND); mpfr_mul(scale[i].v, scale[i].v, tol.v, RND); }
            // Extrapolation tableau, row by row (only the previous row is kept).
            for (int j = 0; j < k; j++) {
                int nsub = nseq[j];
                mpfr_div_ui(h.v, H.v, (unsigned long)nsub, RND);
                mpfr_mul_ui(h2.v, h.v, 2, RND);
                for (int i = 0; i < N; i++) mpfr_set(zprev[i].v, y[i].v, RND);
                axpy(z, y, h, f0);
                for (int s = 1; s < nsub; s++) {
                    sys.rhs(z, fz); evals++;
                    axpy(znext, zprev, h2, fz);
                    std::swap(zprev, z); std::swap(z, znext);
                }
                sys.rhs(z, fz); evals++;
                // Gragg smoothing: (z_n + z_{n-1} + h f(z_n)) / 2
                for (int i = 0; i < N; i++) {
                    mpfr_fma(tmp.v, h.v, fz[i].v, zprev[i].v, RND);
                    mpfr_add(tmp.v, tmp.v, z[i].v, RND);
                    mpfr_mul(row[0][i].v, tmp.v, half.v, RND);
                }
                for (int l = 1; l <= j; l++)
                    for (int i = 0; i < N; i++) {
                        mpfr_sub(tmp.v, row[l - 1][i].v, prev[l - 1][i].v, RND);
                        mpfr_fma(row[l][i].v, tmp.v, coef[j][l].v, row[l - 1][i].v, RND);
                    }
                if (j < k - 1) for (int l = 0; l <= j; l++) std::swap(row[l], prev[l]);
            }
            // Error estimate of the accepted value from the last two diagonal entries.
            double worst = 0.0;
            for (int i = 0; i < N; i++) {
                mpfr_sub(tmp.v, row[k - 1][i].v, row[k - 2][i].v, RND);
                mpfr_abs(tmp.v, tmp.v, RND);
                mpfr_div(ratio.v, tmp.v, scale[i].v, RND);
                if (!mpfr_number_p(ratio.v)) throw std::runtime_error("nonfinite error estimate / collision candidate");
                double rr = mpfr_get_d(ratio.v, RND);
                if (!std::isfinite(rr)) rr = 1e300;  // finite in MPFR but beyond binary64: reject
                if (rr > worst) worst = rr;
            }
            double fac = worst > 0 ? 0.94 * std::pow(0.65 / worst, expo) : 4.0;
            if (fac > 4.0) fac = 4.0;
            if (fac < 0.2) fac = 0.2;
            if (worst <= 1.0) {
                for (int i = 0; i < N; i++) {
                    mpfr_set(y[i].v, row[k - 1][i].v, RND);
                    if (!mpfr_number_p(y[i].v)) throw std::runtime_error("nonfinite state / collision candidate");
                }
                if (last) mpfr_set(t.v, T.v, RND); else mpfr_add(t.v, t.v, H.v, RND);
                steps++;
                if (steps % stride == 0 || mpfr_cmp(t.v, T.v) >= 0) emit(t, y);
                if (!last) { mpfr_mul_d(H.v, H.v, fac, RND); if (mpfr_cmp(H.v, Hmax.v) > 0) mpfr_set(H.v, Hmax.v, RND); }
            } else {
                rejected++;
                mpfr_mul_d(H.v, H.v, fac, RND);
            }
            if (mpfr_cmp_d(H.v, 1e-40) < 0) throw std::runtime_error("step size underflow / collision candidate");
        }
        std::fprintf(out, "# success=1 steps=%ld rejected=%ld rhs_evaluations=%ld\n", steps, rejected, evals);
        std::fclose(out);
        std::printf("success steps=%ld rejected=%ld rhs_evaluations=%ld\n", steps, rejected, evals);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}
