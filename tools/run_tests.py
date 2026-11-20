"""Run the same headless regression checks locally and in GitHub Actions."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import urllib.error

from build_release import ROOT, validate

DEPENDENCIES = {
    "groovy-4.0.28.jar": "org/apache/groovy/groovy/4.0.28/groovy-4.0.28.jar",
    "groovy-json-4.0.28.jar": "org/apache/groovy/groovy-json/4.0.28/groovy-json-4.0.28.jar",
    "ij-1.54p.jar": "net/imagej/ij/1.54p/ij-1.54p.jar",
}


def main():
    validate()
    dependencies = Path(os.environ.get("FIJI_JARS", ROOT / "build/deps"))
    if not os.environ.get("FIJI_JARS"):
        dependencies.mkdir(parents=True, exist_ok=True)
        for name, remote in DEPENDENCIES.items():
            path = dependencies / name
            if not path.is_file():
                temporary = path.with_suffix(".tmp")
                mirrors = ["https://repo.maven.apache.org/maven2/", "https://repo1.maven.org/maven2/", "https://repo.maven.apache.org/maven2/"]
                for attempt in range(3):
                    try:
                        print("Downloading " + mirrors[attempt] + remote, flush=True)
                        with urllib.request.urlopen(mirrors[attempt] + remote, timeout=60) as response:
                            temporary.write_bytes(response.read())
                        break
                    except (urllib.error.URLError, TimeoutError):
                        if attempt == 2:
                            raise
                        time.sleep(2)
                temporary.replace(path)
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
    if "--gui-only" not in sys.argv:
        for test in ["tests/test_progress.groovy", "tests/test_detection_start.groovy", "tests/test_updates.groovy", "tests/test_install.groovy", "tests/test_latest_base.groovy", "tests/test_analysis_bridge.groovy", "tests/test_measurement_protocol.groovy"]:
            subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", test], cwd=ROOT, env=env, check=True)
        subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"], cwd=ROOT, env=env, check=True)
        subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", "tests/test_bootstrap.groovy"], cwd=ROOT, env=env, check=True)
    if "--gui" in sys.argv or "--gui-only" in sys.argv:
        for test in ["tests/test_gui.groovy", "tests/test_layout.groovy", "tests/test_cell_editing.groovy", "tests/test_analysis_gui.groovy"]:
            subprocess.run([java, "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", test], cwd=ROOT, env=env, check=True, timeout=60)


if __name__ == "__main__":
    main()
