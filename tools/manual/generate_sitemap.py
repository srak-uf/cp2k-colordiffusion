#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
from os import path
from datetime import datetime

# Generates a sitemap to guide for search engine crawlers
# See http://www.sitemaps.org for more information.
#
# author: Ole Schuett

def main():
    rootdir = "/var/www/cp2k.org/manual/"

    entries = list()
    for dirname, subdirs, files in os.walk(rootdir, followlinks=True):
        reldir = dirname[len(rootdir):]
        changefreq = "never" if("branch" in reldir.split("/")[0]) else "daily"
        #print('Found directory: %s' % reldir)
        for fn in files:
            if(not fn.lower().endswith(".html")):
                continue # skip non-html files
            absfn = path.join(dirname, fn)
            t = path.getmtime(absfn)
            d = datetime.utcfromtimestamp(t)
            lastmod = d.strftime("%Y-%m-%dT%H:%M:%SZ")
            relfn = path.join(reldir, fn)
            entry  = "<url>"
            entry += "<loc>https://manual.cp2k.org/%s</loc>"%relfn
            entry += "<lastmod>%s</lastmod>"%lastmod
            entry += "<changefreq>%s</changefreq>"%changefreq
            entry += "</url>"
            entries.append(entry)

    print("Found %d html-files."%len(entries))

    fn_out = "sitemap.xml.gz"
    f = gzip.open("sitemap.xml.gz", "wb")
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n')
    f.write("\n".join(entries))
    f.write('\n</urlset>\n')
    f.close()
    print("Wrote "+fn_out)

#===============================================================================
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
    pass #TODO implement selftest
else:
    main()
#EOF
