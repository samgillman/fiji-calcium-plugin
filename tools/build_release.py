"""Build the installable release; validate versions before any publishing."""
from pathlib import Path
import hashlib
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = Path("scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy")


def validate():
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"\d{4}\.\d{1,2}\.\d{1,2}\.\d{1,4}", version):
        raise ValueError("VERSION must have the form YYYY.MM.DD.N")
    source = (ROOT / PLUGIN).read_text(encoding="utf-8")
    if f'static final String VERSION = "{version}"' not in source:
        raise ValueError("CalciumRelease.VERSION does not match VERSION")
    displayed = re.findall(r"version (\d[\d.-]+)", source)
    if not displayed or any(value != version for value in displayed):
        raise ValueError("A plugin version header, title, or log does not match VERSION")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    section = re.search(r"^## " + re.escape(version) + r"\s*\n(.*?)(?=^## |\Z)", changelog, re.M | re.S)
    if not section or not section.group(1).strip():
        raise ValueError("Add release notes for this VERSION to CHANGELOG.md")
    return version, section.group(1).strip()


def build():
    version, notes = validate()
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    source = (ROOT / PLUGIN).read_text(encoding="utf-8")
    shared = source[source.index("class CalciumRelease {"):source.index("/** Background requests")]
    installer = dist / "Install_Calcium_Plugin.groovy"
    installer.write_text("// Fiji Calcium Plugin installer: open in Fiji and press Run.\n" + shared +
                         (ROOT / "tools/installer_entry.groovy").read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    contents = [PLUGIN, Path("README.md"), Path("VERSION"), Path("CHANGELOG.md"), Path("docs/USER_GUIDE.txt")]
    archive = dist / "Calcium_Imaging_Fiji_Plugin.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipped:
        for relative in contents:
            # Fixed metadata makes repeated builds of the same files byte-identical.
            info = zipfile.ZipInfo(relative.as_posix(), date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            zipped.writestr(info, (ROOT / relative).read_bytes())
    plugin_copy = dist / PLUGIN.name
    plugin_copy.write_bytes((ROOT / PLUGIN).read_bytes())
    version_copy = dist / "VERSION.txt"
    version_copy.write_text(version + "\n", encoding="utf-8", newline="\n")
    checksums = "".join(hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.name + "\n" for p in [archive, plugin_copy, version_copy, installer])
    (dist / "SHA256SUMS.txt").write_text(checksums, encoding="utf-8", newline="\n")
    (dist / "release-notes.md").write_text(notes + "\n\n**Install without copying files:** download **Install_Calcium_Plugin.groovy**, open it in Fiji's Script Editor, and press **Run**. It downloads, verifies, backs up and installs the latest plugin. Save your work and restart Fiji. Future updates use **Install update** inside the plugin. The ZIP remains available for manual installation.\n\nChecks use simulated Cellpose runs; this is not a claim of full Windows imaging validation.\n", encoding="utf-8", newline="\n")
    print(f"Built version {version}: {archive}")


if __name__ == "__main__":
    build()
