#!/usr/bin/env python3
"""Verify local guide anchors, table semantics and provenance reading routes."""
from html.parser import HTMLParser
from pathlib import Path
import argparse
class Page(HTMLParser):
 def __init__(self,text):
  super().__init__(); self.ids=set();self.links=[];self.regions=0;self.headers=0;self.feed(text)
 def handle_starttag(self,tag,attrs):
  a=dict(attrs)
  if 'id' in a:self.ids.add(a['id'])
  if tag=='a' and 'href' in a:self.links.append(a['href'])
  if 'linf-table' in a.get('class','').split():
   assert a.get('role')=='region' and a.get('tabindex')=='0' and a.get('aria-label');self.regions+=1
  if tag=='th':self.headers+=1
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args()
for name in ['function-guide','example-datasets']:
 path=a.directory/f'{name}.html';text=path.read_text();page=Page(text)
 assert page.regions>=2 and page.headers>=6,(name,page.regions,page.headers)
 for link in page.links:
  if link.startswith('#') and len(link)>1: assert link[1:] in page.ids,(name,link)
 assert 'table-layout: fixed' not in text and 'body { max-width: 960px' not in text
 print(f'{name}: local anchors and {page.regions} named keyboard-focusable tables verified')
