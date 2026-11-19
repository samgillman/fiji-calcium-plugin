"""Bundle the exact pinned R sources into the single-file plugin, preserving older updaters."""
from pathlib import Path
import base64, hashlib, io, zipfile
ROOT=Path(__file__).resolve().parents[1]
PLUGIN=ROOT/'scripts/Plugins/Calcium Imaging/Calcium_Cells_Step_By_Step.groovy'
START='// BEGIN GENERATED SIMPLECA BUNDLE'
END='// END GENERATED SIMPLECA BUNDLE'
def bundle():
    out=io.BytesIO()
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted((ROOT/'analysis').rglob('*')):
            if p.is_file() and (p.suffix in ['.R','.tsv'] or p.name=='LICENSE'):
                info=zipfile.ZipInfo(p.relative_to(ROOT/'analysis').as_posix(),(2026,1,1,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED
                z.writestr(info,p.read_bytes())
    data=out.getvalue()
    return START+'\nclass CalciumSimpleCaBundle {\n    static final String HASH = "'+hashlib.sha256(data).hexdigest()+'"\n    static final String DATA = "'+base64.b64encode(data).decode()+'"\n}\n'+END

def sync(check=False):
    s=PLUGIN.read_text(encoding='utf-8'); expected=bundle()
    if START not in s: raise ValueError('Missing bundle marker')
    a=s.index(START); b=s.index(END,a)+len(END)
    if check:
        if s[a:b]!=expected: raise ValueError('Run python tools/embed_analysis.py before packaging')
    else: PLUGIN.write_text(s[:a]+expected+s[b:],encoding='utf-8',newline='\n')
if __name__=='__main__': sync()
