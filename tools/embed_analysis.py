"""Bundle the exact pinned R sources into the single-file plugin, preserving older updaters."""
from pathlib import Path
import base64, hashlib, io, re, zipfile
ROOT=Path(__file__).resolve().parents[1]
PLUGIN=ROOT/'scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'
START='// BEGIN GENERATED SIMPLECA BUNDLE'
END='// END GENERATED SIMPLECA BUNDLE'
def sources():
    return {p.relative_to(ROOT/'analysis').as_posix():p.read_bytes()
            for p in sorted((ROOT/'analysis').rglob('*'))
            if p.is_file() and (p.suffix in ['.R','.tsv'] or p.name=='LICENSE')}
def bundle():
    out=io.BytesIO()
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        for name,content in sources().items():
            info=zipfile.ZipInfo(name,(2026,1,1,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED
            info.create_system=3; info.external_attr=0o100644 << 16
            z.writestr(info,content)
    data=out.getvalue()
    return START+'\nclass CalciumSimpleCaBundle {\n    static final String HASH = "'+hashlib.sha256(data).hexdigest()+'"\n    static final String DATA = "'+base64.b64encode(data).decode()+'"\n}\n'+END

# The Step 4 web page travels inside the same single plugin file, in a block of its own.
WEB_START='// BEGIN GENERATED WEB BUNDLE'
WEB_END='// END GENERATED WEB BUNDLE'
def web_sources():
    return {p.name:p.read_bytes() for p in sorted((ROOT/'web').glob('step4.*')) if p.suffix in ['.html','.css','.js']}
def web_bundle():
    out=io.BytesIO()
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        for name,content in web_sources().items():
            info=zipfile.ZipInfo(name,(2026,1,1,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED
            info.create_system=3; info.external_attr=0o100644 << 16
            z.writestr(info,content)
    data=out.getvalue()
    return WEB_START+'\nclass CalciumWebBundle {\n    static final String SHA256 = "'+hashlib.sha256(data).hexdigest()+'"\n    static final String CONTENT = "'+base64.b64encode(data).decode()+'"\n}\n'+WEB_END
def sync_web(s,check):
    if WEB_START not in s: raise ValueError('Missing web bundle marker')
    a=s.index(WEB_START); b=s.index(WEB_END,a)+len(WEB_END)
    if not check: return s[:a]+web_bundle()+s[b:]
    digest=re.search(r'SHA256 = "([a-f0-9]{64})"',s[a:b]).group(1)
    data=base64.b64decode(re.search(r'CONTENT = "([A-Za-z0-9+/=]+)"',s[a:b]).group(1),validate=True)
    if hashlib.sha256(data).hexdigest()!=digest: raise ValueError('Invalid embedded Step 4 page checksum')
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        expected=web_sources()
        if len(z.namelist())!=len(expected) or {name:z.read(name) for name in z.namelist()}!=expected:
            raise ValueError('Run python tools/embed_analysis.py before packaging')
    return s

def sync(check=False):
    s=PLUGIN.read_text(encoding='utf-8')
    if START not in s: raise ValueError('Missing bundle marker')
    a=s.index(START); b=s.index(END,a)+len(END)
    if check:
        # ZIP host metadata and zlib output can differ across OS/Python builds.
        # Verify the installed payload's digest and exact uncompressed sources instead.
        digest=re.search(r'HASH = "([a-f0-9]{64})"',s[a:b]).group(1)
        data=base64.b64decode(re.search(r'DATA = "([A-Za-z0-9+/=]+)"',s[a:b]).group(1),validate=True)
        if hashlib.sha256(data).hexdigest()!=digest: raise ValueError('Invalid embedded analysis checksum')
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            expected=sources()
            if len(z.namelist())!=len(expected) or {name:z.read(name) for name in z.namelist()}!=expected:
                raise ValueError('Run python tools/embed_analysis.py before packaging')
        sync_web(s,True)
    else: PLUGIN.write_text(sync_web(s[:a]+bundle()+s[b:],False),encoding='utf-8',newline='\n')
if __name__=='__main__': sync()
