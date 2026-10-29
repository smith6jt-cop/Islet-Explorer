/* Read-only verification after derived context/parenting refresh. */
import com.google.gson.Gson
def cfg = new Gson().fromJson(new File(args[0]).getText('UTF-8'), Map)
def invariant = evaluate(new File(cfg.invariants_script))
def checks = invariant.verify(getCurrentHierarchy(), invariant.readJSON(new File(cfg.inventory)))
def ids = invariant.readJSON(new File(cfg.new_ids)) as Set
def current = new HashSet(getDetectionObjects().collect { it.getID().toString() })
def expected = new HashSet(invariant.readJSON(new File(cfg.inventory)).cells.keySet())
expected.addAll(ids)
// Groovy's ==/!= coerces Set elements with a quadratic pairwise comparison. UUID strings
// need ordinary hash membership; keep exact equality practical for hundreds of thousands of cells.
if (current.size() != expected.size() || !current.containsAll(expected))
    throw new IllegalStateException('Cell inventory differs from original plus accepted additions')
V3_RESULT = [ok: true, changed: false, message: 'Original objects preserved and all new IDs present', details: checks]
