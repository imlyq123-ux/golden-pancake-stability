"""Newton-correct the supplied mu=0.613 seed with the Taylor/Decimal pipeline (85 digits) and compare with GBS."""
import json, sys
from decimal import getcontext
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import high_precision as hp
getcontext().prec = 85
out = Path(__file__).resolve().parent
p, hist = hp.correct('0.613', hp.SEED, 85, 50, '1e-64', out, 'taylor85_0p613', target='1e-66', stride=10**9)
dat = hp.verify_orbit('0.613', p, hist, 85, 56, '1e-68', out, 'taylor85_0p613', stride=20)
(out / 'taylor85_0p613.json').write_text(json.dumps(dat, indent=1))
print('closure', dat['max_closure'], 'u', dat['u'], 'v', dat['v'], 'T', dat['T'])
