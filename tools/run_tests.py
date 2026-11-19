"""Run the same headless regression checks locally and in GitHub Actions."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import urllib.request

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
                with urllib.request.urlopen("https://repo.maven.apache.org/maven2/" + remote, timeout=60) as response:
                    temporary.write_bytes(response.read())
                temporary.replace(path)
    jars = [dependencies / name for name in DEPENDENCIES]
    for jar in jars:
        if not jar.is_file():
            raise FileNotFoundError(jar)
    java = str(Path(os.environ["JAVA_HOME"]) / "bin" / ("java.exe" if os.name == "nt" else "java")) if os.environ.get("JAVA_HOME") else shutil.which("java")
    if not java:
        raise RuntimeError("Java 17 or newer is required")
    env = dict(os.environ, CALCIUM_TEST_PYTHON=sys.executable)
    for test in ["tests/test_progress.groovy", "tests/test_updates.groovy", "tests/test_latest_base.groovy"]:
        subprocess.run([java, "-Djava.awt.headless=true", "-Dfile.encoding=UTF-8", "-cp", os.pathsep.join(str(jar) for jar in jars), "groovy.ui.GroovyMain", test], cwd=ROOT, env=env, check=True)
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"], cwd=ROOT, env=env, check=True)


if __name__ == "__main__":
    main()
