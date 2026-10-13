"""A full-period, independently initialized time-reversal monodromy check."""
import json
from pathlib import Path
from decimal import Decimal as D, getcontext
from high_precision import verify_orbit, inverse, sub, maxabs, ROOT


def main():
    source=json.loads((ROOT/'results/high_precision_orbits.json').read_text())
    getcontext().prec=source['precision_digits']
    orbit=next(r for r in source['orbits'] if D(r['mu'])==D('.613'))
    p=[-D(orbit['u']),-D(orbit['v']),D(orbit['T'])]
    reversed_orbit=verify_orbit('.613',p,[],65,44,'1e-46',ROOT/'results/mp_raw','time_reversed_seed')
    forward=[[D(x) for x in row] for row in orbit['floquet']['full_monodromy']]
    inv=inverse(forward)
    signs=[D(1)]*6+[D(-1)]*6
    expected=[[signs[i]*inv[i][j]*signs[j] for j in range(12)] for i in range(12)]
    actual=[[D(x) for x in row] for row in reversed_orbit['floquet']['full_monodromy']]
    result={'relation':'Phi_reversed = R * inverse(Phi_forward) * R, R=diag(I6,-I6), up to actual periodic closure error.',
        'monodromy_max_identity_defect':str(maxabs(sub(actual,expected))),
        'forward_source':'results/high_precision_orbits.json','reversed_orbit':reversed_orbit,
        'limitation':'Independent initial velocity signs and full-period integration; same Taylor implementation.'}
    (ROOT/'results/time_reversal_verification.json').write_text(json.dumps(result,indent=2))
    print('Time-reversal monodromy identity defect:',result['monodromy_max_identity_defect'],flush=True)
    assert D(result['monodromy_max_identity_defect'])<D('1e-28')


if __name__=='__main__':main()
