PY ?= python3
.PHONY: verify verify-full
verify:        ## standard library only
	$(PY) -I verify.py CLAIMS.md
verify-full:   ## also regenerates every number and table of the paper (needs numpy + scipy: pip install -r requirements.txt)
	$(PY) -I verify.py CLAIMS.md --full
