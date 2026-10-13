"""Half-period reflection symmetry of the Golden Pancake orbits and the half-period map N.

Each Euler-section orbit on the isola satisfies X(T/2) = sigma X(0), where sigma reflects every
position and velocity across a line through the centre of mass (axis angle depends on mu).
Reflections commute with the Newtonian flow, hence X(t + T/2) = sigma X(t) and the monodromy is
a square, M = N^2 with N = sigma * DPhi_{T/2}(X(0)). Consequences used in the analysis:
  * nontrivial pair sums satisfy s_M = s_N^2 - 2, so s_M + 2 = s_N^2 >= 0 when s_N is real:
    M-multipliers can only touch -1 (when N-multipliers pass +-i) -- no period doubling of M
    inside the symmetric family; and |lambda_M| = |lambda_N|^2, so M and N are spectrally
    stable together;
  * N has trivial multipliers +1 (time/energy, translation along the axis) and -1
    (rotation/angular momentum, translation across the axis), four each.
This script integrates X(0) -> X(T/2) with the 12 variational columns using the independent
GBS integrator, measures the symmetry defect, forms N, compares N^2 with the saved full-period
monodromy, and reports the nontrivial spectrum of N.
"""
import argparse
import json
from pathlib import Path

import mpmath
from mpmath import mp, mpf

import independent_mp as im

ROOT = Path(__file__).resolve().parent


def analyse(record, digits, tol, stages, raw_dir):
    mu = record['mu']
    T = mpf(record['T'])
    half = im.s(T / 2)
    path = raw_dir / f"half_{mu.replace('.', 'p')[:60]}.csv"
    rows, info = im.run_gbs(mu, [record['u'], record['v'], half], digits, tol, stages, 12, path, 10 ** 9)
    final = im.parse(rows[-1])
    x0 = im.initial_state(mu, record['u'], record['v'])
    xh = final[1:13]
    c2, s2 = -xh[0], -xh[1]                     # sigma(-1, 0) = (-cos 2a, -sin 2a) = q1(T/2)
    norm = mpmath.sqrt(c2 ** 2 + s2 ** 2)
    c2, s2 = c2 / norm, s2 / norm
    sigma = mpmath.zeros(12, 12)
    for b in range(6):
        sigma[2 * b, 2 * b], sigma[2 * b, 2 * b + 1] = c2, s2
        sigma[2 * b + 1, 2 * b], sigma[2 * b + 1, 2 * b + 1] = s2, -c2
    sx0 = sigma * mpmath.matrix(x0)
    defect = max(abs(sx0[i] - xh[i]) for i in range(12))
    A = mpmath.matrix(12, 12)
    for c in range(12):
        for r in range(12):
            A[r, c] = final[13 + 12 * c + r]
    N = sigma * A
    M = mpmath.matrix([[mpf(x) for x in row] for row in record['full_monodromy']])
    square_defect = max(abs(x) for x in (N * N - M))
    traces, power = [], N
    for k in range(1, 5):
        traces.append(sum(power[i, i] for i in range(12)))
        power = power * N
    p = [traces[k - 1] - 4 - 4 * (-1) ** k for k in range(1, 5)]
    e1 = p[0]
    e2 = (e1 * p[0] - p[1]) / 2
    disc = 2 * p[1] + 8 - p[0] ** 2
    root = mpmath.sqrt(mpmath.mpc(disc))
    s_n = [(e1 + root) / 2, (e1 - root) / 2]
    s_m = [mpf(x[0]) for x in record['floquet']['pair_sums_s']]
    implied = sorted(mpmath.re(z * z - 2) for z in s_n)
    return {'mu': mu, 'axis_angle_degrees': mpmath.nstr(mpmath.degrees(mpmath.atan2(s2, c2) / 2) % 180, 20),
            'symmetry_defect_max_abs': mpmath.nstr(defect, 6), 'N_squared_minus_M_max_abs': mpmath.nstr(square_defect, 6),
            'D_N': mpmath.nstr(disc, 25), 's_N': [[mpmath.nstr(mpmath.re(z), 25), mpmath.nstr(mpmath.im(z), 25)] for z in s_n],
            's_M_saved': [mpmath.nstr(x, 25) for x in sorted(s_m)],
            's_N_squared_minus_2': [mpmath.nstr(x, 25) for x in implied],
            'pair_sum_identity_defect': mpmath.nstr(max(abs(a - b) for a, b in zip(sorted(s_m), implied)), 6),
            'half_period_integration': info}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--input', type=Path, default=ROOT / 'results' / 'independent_gbs_orbits.json')
    parser.add_argument('--digits', type=int, default=80)
    parser.add_argument('--tol', default='1e-66')
    parser.add_argument('--stages', type=int, default=24)
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'half_period_symmetry.json')
    parser.add_argument('--raw-dir', type=Path, default=ROOT / 'results' / 'gbs_raw_half_period')
    args = parser.parse_args()
    mp.dps = args.digits
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    records = json.loads(args.input.read_text())['orbits']
    out = []
    for rec in records:
        res = analyse(rec, args.digits, args.tol, args.stages, args.raw_dir)
        out.append(res)
        print(f"mu={res['mu'][:22]:22s} axis={res['axis_angle_degrees'][:10]} sym_defect={res['symmetry_defect_max_abs']} "
              f"|N^2-M|={res['N_squared_minus_M_max_abs']} D_N={mpmath.nstr(mpf(res['D_N']), 8)} "
              f"identity_defect={res['pair_sum_identity_defect']}", flush=True)
    args.output.write_text(json.dumps({'description': __doc__, 'source': str(args.input), 'digits': args.digits,
                                       'results': out}, indent=1))


if __name__ == '__main__':
    main()
