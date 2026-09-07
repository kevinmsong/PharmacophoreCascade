"""Render DOCX with native Microsoft Word on Windows, then emit page PNGs for QA.

Word is the installed renderer on this workstation; LibreOffice is unavailable.
Uses a separate hidden Word instance and closes documents without changing them.
"""
import argparse
from pathlib import Path
import fitz
import pythoncom
import win32com.client


def render(paths,out,dpi):
    pythoncom.CoInitialize()
    word=win32com.client.DispatchEx('Word.Application')
    word.Visible=False
    word.DisplayAlerts=0
    try:
        for path in paths:
            folder=out/path.stem;folder.mkdir(parents=True,exist_ok=True)
            pdf=folder/(path.stem+'.pdf')
            document=word.Documents.Open(str(path.resolve()),ReadOnly=True,AddToRecentFiles=False)
            try:
                document.Repaginate()
                document.ExportAsFixedFormat(str(pdf.resolve()),17,OpenAfterExport=False)
            finally:
                document.Close(SaveChanges=0)
            with fitz.open(pdf) as doc:
                for i,page in enumerate(doc,1):
                    page.get_pixmap(dpi=dpi,alpha=False).save(folder/f'page-{i:03d}.png')
                print(f'{path.name}: {len(doc)} rendered pages in {folder}',flush=True)
    finally:
        word.Quit(SaveChanges=0)
        pythoncom.CoUninitialize()


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('documents',nargs='+',type=Path)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--dpi',type=int,default=110)
    args=ap.parse_args()
    render(args.documents,args.output_dir,args.dpi)
