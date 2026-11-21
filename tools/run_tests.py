"""Run the same headless regression checks locally and in GitHub Actions."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import time
import hashlib
import urllib.request
import urllib.error

from build_release import ROOT, validate

DEPENDENCIES = {
    "groovy-4.0.28.jar": "org/apache/groovy/groovy/4.0.28/groovy-4.0.28.jar",
    "groovy-json-4.0.28.jar": "org/apache/groovy/groovy-json/4.0.28/groovy-json-4.0.28.jar",
    "ij-1.54p.jar": "net/imagej/ij/1.54p/ij-1.54p.jar",
}
# SHA-256 of the exact jars on Maven Central, so a mirror or a damaged
# download can never stand in for the pinned version.
CHECKSUMS = {
    "groovy-4.0.28.jar": "d77164a44ab994837e61eecead64721bdb1794fa2bc6704ac6fe70cd32ef9dd3",
    "groovy-json-4.0.28.jar": "edec45c85ff811f9a9baacdf5b020e34593a0055ac63c92fbf1233c3b2cc9cc6",
    "ij-1.54p.jar": "2e1a09961dfb41cee66ddc821b2577a41a072566ce45a49bae69267099741e20",
}
# Maven Central intermittently answers GitHub Actions runners with 429 or a
# transient 404, and a push to main starts two Checks runs at once. Spread
# attempts over two CDNs and a Google-hosted copy of Central, backing off.
MIRRORS = [
    "https://repo.maven.apache.org/maven2/",
    "https://repo1.maven.org/maven2/",
    "https://maven-central.storage-download.googleapis.com/maven2/",
]
ATTEMPTS = 8


def download(name, remote, path):
    temporary = path.with_suffix(".tmp")
    last_error = None
    for attempt in range(ATTEMPTS):
        url = MIRRORS[attempt % len(MIRRORS)] + remote
        try:
            print("Downloading " + url, flush=True)
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            digest = hashlib.sha256(data).hexdigest()
            if digest != CHECKSUMS[name]:
                raise ValueError(name + " from " + url + " has SHA-256 " + digest + ", expected " + CHECKSUMS[name])
            temporary.write_bytes(data)
            temporary.replace(path)
            return
        except (urllib.error.URLError, TimeoutError, ConnectionError, ValueError) as error:
            last_error = error
            print("  failed: " + str(error), flush=True)
            if attempt == ATTEMPTS - 1:
                break
            delay = min(2 ** attempt, 30)
            retry_after = getattr(error, "headers", None) and error.headers.get("Retry-After")
            if retry_after and retry_after.isdigit():
                delay = min(max(delay, int(retry_after)), 60)
            time.sleep(delay)
    raise RuntimeError("Could not download " + name + " after " + str(ATTEMPTS) + " attempts") from last_error


def main():
    validate()
    dependencies = Path(os.environ.get("FIJI_JARS", ROOT / "build/deps"))
    if not os.environ.get("FIJI_JARS"):
        dependencies.mkdir(parents=True, exist_ok=True)
        for name, remote in DEPENDENCIES.items():
            path = dependencies / name
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != CHECKSUMS[name]:
                download(name, remote, path)
    jars = [dependencies / name for name in DEPENDENCIES]
    for jar in jars:
        if not jar.is_file():
            raise FileNotFoundError(jar)
    java = str(Path(os.environ["JAVA_HOME"]) / "bin" / ("java.exe" if os.name == "nt" else "java")) if os.environ.get("JAVA_HOME") else shutil.which("java")
    if not java:
        raise RuntimeError("Java 17 or newer is required")
    env = dict(os.environ, CALCIUM_TEST_PYTHON=sys.executable)
    if "--analysis-portable" in sys.argv:
        subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", "tests/test_analysis_portable.groovy"], cwd=ROOT, env=env, check=True, timeout=300)
        return
    if "--analysis-only" in sys.argv:
        subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", "tests/test_analysis_runtime.groovy"], cwd=ROOT, env=env, check=True, timeout=1200)
        return
    if "--browser" in sys.argv:
        # Drive the real local Step 4 page. Install Playwright and its Chromium
        # runtime first; this gate must not silently fall back to syntax checks.
        node = shutil.which("node")
        if not node:
            raise RuntimeError("The browser check needs Node with Playwright")
        try:
            subprocess.run([node, "-e", "require.resolve('playwright')"], cwd=ROOT, env=env, check=True)
        except subprocess.CalledProcessError as error:
            raise RuntimeError("The browser check needs the Playwright Node package and its Chromium runtime; see docs/RELEASING.md") from error
        folder = ROOT / "build/web-step4"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "url.txt").unlink(missing_ok=True)
        # Settings remembered from an earlier run would change what the page opens with; every run starts as a clean checkout does.
        (folder / "step4_last_settings.json").unlink(missing_ok=True)
        served = subprocess.Popen([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", "tests/browser/serve_step4.groovy"], cwd=ROOT, env=env)
        try:
            deadline = time.time() + 60
            while not (folder / "url.txt").is_file():
                if served.poll() is not None or time.time() > deadline:
                    raise RuntimeError("The Step 4 page did not start")
                time.sleep(0.2)
            subprocess.run([node, "tests/browser/step4_page.js", str(folder)], cwd=ROOT, env=env, check=True, timeout=120)
        finally:
            (folder / "stop").write_text("stop", encoding="utf-8")
            try:
                served.wait(timeout=10)
            except subprocess.TimeoutExpired:
                served.kill()
                served.wait(timeout=10)
        return
    if "--gui-only" not in sys.argv:
        node = shutil.which("node")
        if node:  # the page's script must at least parse; the full browser check is --browser
            subprocess.run([node, "--check", "web/step4.js"], cwd=ROOT, env=env, check=True)
        for test in ["tests/test_progress.groovy", "tests/test_detection_start.groovy", "tests/test_updates.groovy", "tests/test_install.groovy", "tests/test_latest_base.groovy", "tests/test_analysis_bridge.groovy", "tests/test_measurement_protocol.groovy", "tests/test_web_step4.groovy", "tests/test_sensitivity.groovy", "tests/test_step4_carryover.groovy", "tests/test_motion.groovy", "tests/test_grading.groovy"]:
            subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", test], cwd=ROOT, env=env, check=True)
        subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"], cwd=ROOT, env=env, check=True)
        subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", "tests/test_bootstrap.groovy"], cwd=ROOT, env=env, check=True)
    if "--gui" in sys.argv or "--gui-only" in sys.argv:
        for test in ["tests/test_gui.groovy", "tests/test_start_updates.groovy", "tests/test_layout.groovy", "tests/test_cell_editing.groovy", "tests/test_picture_views.groovy", "tests/test_analysis_gui.groovy"]:
            subprocess.run([java, "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", test], cwd=ROOT, env=env, check=True, timeout=60)


if __name__ == "__main__":
    main()
