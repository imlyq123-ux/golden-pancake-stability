"""Compatibility entry point for the portable boundary CLI."""
import sys
from pancake_stability_boundaries import main
if __name__=='__main__':
    sys.argv.append('--upper-only')
    main()
