#!/usr/bin/env python3

from dataclasses import dataclass, field
from typing import Optional, Any
from pathlib import Path
from shutil import copytree
from sys import argv
import inspect
import json
import jinja2
import vipster

if len(argv) > 1:
    targetdir = Path(argv[1])
    if targetdir.exists():
        if not targetdir.is_dir():
            raise ValueError("Path {} exists but is not a directory".format(targetdir))
    else:
        targetdir.mkdir()
else:
    targetdir = Path('.')

#########################################
# Extract Type hierarchy and docstrings #
#########################################

def getTypeStr(t, parent):
    if inspect.ismodule(t):
        return "Module"
    if inspect.isclass(t):
        return "Class"
    if inspect.isroutine(t):
        return "Function"
    if type(t) == property:
        return "Property"
    if type(t) == parent:
        return "Instance"
    return str(type(t))


@dataclass
class Node:
    obj: Any
    type: str
    name: str
    doc: str
    parent: Optional['Node']
    children: ['Node'] = field(default_factory=list)

root = Node(vipster, "Module", "vipster", "", None)

def walkTree(node):
    # TODO: get members of Enums and _SOME_ instances?
    if node.type not in ["Module", "Class"]:
        return
    for name, child in inspect.getmembers(node.obj):
        if name[0] == '_':
            continue
        node.children.append(Node(child,
                                  getTypeStr(child, node.obj),
                                  name,
                                  doc if (doc := inspect.getdoc(child)) else "",
                                  node))
        walkTree(node.children[-1])

walkTree(root)

########################
# Extract Format plugins
########################

formats = [c for n,c in inspect.getmembers(vipster.Plugins) if n[0] != '_']

###############
# Create Page #
###############

copytree(str(Path(__file__).parent.absolute())+"/gh-pages",
         targetdir, dirs_exist_ok=True)

env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(Path(__file__).parent.absolute()),
    autoescape=jinja2.select_autoescape()
)

@dataclass
class Page:
    name: str
    url: str
    base: str = ""

about_content = ""
qt_content = ""
download_content = ""

all_pages = [Page("About", "about.html"),
             Page("Documentation", "gui.html"),
             Page("GUI Features", "gui.html", "Documentation"),
             Page("File Formats", "format.html", "Documentation"),
             Page("Python API", "python.html", "Documentation"),
             Page("Download", "download.html"),
            ]

def printPages(pages):
    for page in pages:
        template = env.get_template(page.url)
        with open(str(targetdir.absolute())+"/"+page.url, 'w') as f:
            f.write(template.render(current=page, pages=all_pages, tree=[root], formats=formats))

printPages(all_pages)
