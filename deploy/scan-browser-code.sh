#!/usr/bin/env bash
# Local OSS SAST only. Fetching the public rule snapshot never sends source.
set -euo pipefail

usage() { echo "Usage: $0 --output NEW_DIRECTORY [--rules LOCAL_YAML]" >&2; exit 2; }
output=''
rules=''
while (( $# )); do
  case "$1" in
    --output) (( $# >= 2 )) || usage; output=$2; shift 2 ;;
    --rules) (( $# >= 2 )) || usage; rules=$2; shift 2 ;;
    *) usage ;;
  esac
done
[[ -n "$output" && ! -e "$output" ]] || usage
repo_root=$(git rev-parse --show-toplevel)
semgrep_bin=${SEMGREP_BIN:-semgrep}
for dependency in Rscript node python3 curl sha256sum; do command -v "$dependency" >/dev/null; done
command -v "$semgrep_bin" >/dev/null
mkdir -p "$output"
output=$(cd "$output" && pwd)
export SEMGREP_SETTINGS_FILE="$output/semgrep-settings.yml"
export SEMGREP_LOG_FILE="$output/semgrep-internal.log"
export SEMGREP_ENABLE_VERSION_CHECK=0 SEMGREP_SEND_METRICS=off
unset SEMGREP_APP_TOKEN SEMGREP_RULES
"$semgrep_bin" --disable-version-check --version > "$output/semgrep-version.txt" 2> "$output/semgrep-version-stderr.txt"
[[ $(cat "$output/semgrep-version.txt") == '1.176.1' ]] || { echo 'Semgrep 1.176.1 is required.' >&2; exit 1; }
rule_url='https://semgrep.dev/c/p/javascript'
if [[ -n "$rules" ]]; then
  rule_acquisition='local_snapshot_replay'
  cp "$rules" "$output/javascript-rules.yml"
else
  rule_acquisition='official_registry_download'
  curl --fail --location --proto '=https' --tlsv1.2 --silent --show-error "$rule_url" --output "$output/javascript-rules.yml"
fi
sha256sum "$output/javascript-rules.yml" | awk '{print $1}' > "$output/rules.sha256"
mkdir -p "$output/targets"
Rscript "$repo_root/scripts/extract_browser_scripts.R" "$repo_root/app/shiny_app" "$output/targets/embedded-r" > "$output/extraction.log" 2>&1

python3 - "$repo_root" "$output" "$rule_url" "$rule_acquisition" <<'PY'
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from html.parser import HTMLParser

repo, evidence = map(Path, sys.argv[1:3])
www = repo / "app/shiny_app/www"
manifest = []

class InlineScripts(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.active = None
        self.scripts = []
    def handle_starttag(self, tag, attrs):
        if tag != "script":
            return
        attrs = dict(attrs)
        if not attrs.get("src") and attrs.get("type", "").lower() in ("", "module", "text/javascript", "application/javascript"):
            self.active = [self.getpos()[0], attrs.get("type") == "module", ""]
    def handle_data(self, data):
        if self.active is not None:
            self.active[2] += data
    def handle_endtag(self, tag):
        if tag == "script" and self.active is not None:
            if self.active[2].strip():
                self.scripts.append(self.active)
            self.active = None

for source in sorted(www.rglob("*")):
    if not source.is_file() or source.relative_to(www).parts[:2] == ("avivator", "assets"):
        continue
    suffix = source.suffix.lower()
    if suffix not in (".js", ".mjs", ".html"):
        continue
    relative = source.relative_to(repo)
    manifest.append({"origin": str(relative), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
    if suffix in (".js", ".mjs"):
        target = evidence / "targets" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    else:
        parser = InlineScripts()
        parser.feed(source.read_text(encoding="utf-8"))
        if parser.active is not None:
            raise SystemExit("Unclosed inline script element")
        for line, module, code in parser.scripts:
            target = evidence / "targets" / relative.with_suffix(f".L{line}.{'mjs' if module else 'js'}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(code, encoding="utf-8")
for source in sorted((repo / "app/shiny_app/R").glob("*.R")) + [repo / "app/shiny_app/app.R"]:
    manifest.append({"origin": str(source.relative_to(repo)), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
report = {
    "schema_version": 1,
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
    "git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo)),
    "semgrep_version": "1.176.1",
    "rule_set": "p/javascript",
    "rule_source": sys.argv[3],
    "rule_acquisition": sys.argv[4],
    "rule_snapshot_sha256": (evidence / "rules.sha256").read_text().strip(),
    "sources": manifest,
    "scope": "Authored deployed JS/mjs, inline HTML scripts and parsed browser literals in application R.",
    "limitations": ["OSS JavaScript rules are scoped SAST, not a complete security assessment.",
                    "Synthetic R format values cannot establish cross-language taint or runtime safety.",
                    "Dynamically assembled fragments and package-generated scripts are outside extraction coverage.",
                    "Vendored Avivator assets are excluded and require the separate vendor inventory/Retire review."],
    "execution": "Local OSS scan; metrics/version checks/tracing/remote validation disabled; source is not uploaded.",
}
(evidence / "scan-provenance.json").write_text(json.dumps(report, indent=2) + "\n")
PY

mapfile -d '' targets < <(find "$output/targets" -type f \( -name '*.js' -o -name '*.mjs' \) -print0)
(( ${#targets[@]} > 0 )) || { echo 'No browser scan targets found.' >&2; exit 1; }
for target in "${targets[@]}"; do
  node --check "$target" >> "$output/syntax-check.log" 2>&1
done
scan_status=0
"$semgrep_bin" scan --config "$output/javascript-rules.yml" \
  --metrics=off --disable-version-check --no-trace --no-secrets-validation \
  --oss-only --strict --error --disable-nosem --no-git-ignore --max-target-bytes 0 \
  --json --output "$output/semgrep.json" "$output/targets" \
  > "$output/semgrep-stdout.log" 2> "$output/semgrep-stderr.log" || scan_status=$?
python3 - "$output" "$scan_status" "${#targets[@]}" <<'PY'
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
status, expected = map(int, sys.argv[2:])
report = json.loads((root / "semgrep.json").read_text())
scanned = report.get("paths", {}).get("scanned", [])
summary = {"semgrep_exit_status": status, "expected_targets": expected,
           "scanned_targets": len(scanned), "findings": len(report.get("results", [])),
           "errors": len(report.get("errors", [])),
           "skipped_rules": len(report.get("skipped_rules", [])),
           "passed": status == 0 and len(scanned) == expected and not report.get("results") and not report.get("errors") and not report.get("skipped_rules")}
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary))
raise SystemExit(0 if summary["passed"] else 1)
PY
