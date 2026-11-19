"""Decide whether a tested version needs publishing, without overwriting releases."""
import json
import os
import re
import subprocess
from pathlib import Path

from build_release import validate


def parts(version):
    if not re.fullmatch(r"v?\d{4}\.\d{1,2}\.\d{1,2}\.\d{1,4}", version):
        raise ValueError("Unrecognised release version: " + version)
    return tuple(map(int, version.removeprefix("v").split(".")))


def release(repository, suffix, run=subprocess.run):
    result = run(["gh", "api", f"repos/{repository}/releases/{suffix}"],
                 capture_output=True, text=True, check=False)
    try:
        data = json.loads(result.stdout)
    except (ValueError, TypeError) as error:
        raise RuntimeError("Could not read GitHub release metadata: " + result.stderr) from error
    if result.returncode:
        if str(data.get("status")) == "404" and data.get("message") == "Not Found":
            return None
        raise RuntimeError("GitHub release lookup failed: " + result.stderr)
    return data


def needed(version, repository, lookup=release):
    existing = lookup(repository, "tags/v" + version)
    if existing is not None:
        if existing.get("draft") or existing.get("prerelease"):
            raise RuntimeError("This version already has a draft/prerelease; review it before publishing")
        return False
    latest = lookup(repository, "latest")
    if latest and parts(version) <= parts(latest["tag_name"]):
        raise RuntimeError("The new VERSION must be higher than the latest published release")
    return True


def main():
    version, _ = validate()
    if os.environ.get("GITHUB_REF_TYPE") == "tag" and os.environ["GITHUB_REF_NAME"] != "v" + version:
        raise ValueError("Tag must match VERSION")
    publish = needed(version, os.environ["GITHUB_REPOSITORY"])
    with Path(os.environ["GITHUB_OUTPUT"]).open("a", encoding="utf-8") as output:
        output.write("needed=" + str(publish).lower() + "\n")
    print(f"Version {version}: " + ("ready to publish" if publish else "already published; leaving release unchanged"))


if __name__ == "__main__":
    main()
