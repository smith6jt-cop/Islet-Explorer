/* Reuse accepted distance/parenting code, then verify preserved raw objects before the footer saves. */
import com.google.gson.Gson
def augmentationConfig = new Gson().fromJson(new File(args[0]).getText('UTF-8'), Map)
def augmentationHelpers = evaluate(new File(augmentationConfig.invariants_script))
String augmentationMode = args[1]
if (!(augmentationMode in ['3', '8'])) throw new IllegalArgumentException('Expected context step 3 or 8')
String augmentationSource = augmentationMode == '3' ? '3_signed_distance.groovy' : '8_hierarchy_parenting.groovy'
// Nested Groovy evaluate() does not inherit QuPath's editor/compiler static imports.
// Supply the same scripting API explicitly so PROJECT_BASE_DIR and helper methods resolve.
def augmentationSourceFile = new File(new File(augmentationConfig.invariants_script).getParentFile(), augmentationSource)
evaluate('import static qupath.lib.gui.scripting.QPEx.*\n' + augmentationSourceFile.getText('UTF-8'))
def augmentationChecks = augmentationHelpers.verify(getCurrentHierarchy(),
    augmentationHelpers.readJSON(new File(augmentationConfig.inventory)))
augmentationHelpers.requireSaveWindow(augmentationConfig)
V3_RESULT = [ok: true, changed: true, message: "Refreshed ${augmentationMode}; original raw objects preserved",
             details: augmentationChecks]
