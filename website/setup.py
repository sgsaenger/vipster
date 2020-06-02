#!/usr/bin/env python3

from dataclasses import dataclass, field
from typing import Optional, Any
from pathlib import Path
import inspect
import vipster
import json
import jinja2
from sys import argv

if argv[1]:
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

###############
# Create Page #
###############

env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(Path(__file__).parent.absolute()),
    autoescape=jinja2.select_autoescape()
)

@dataclass
class Page:
    name: str
    url: str
    pages: ['Page'] = field(default_factory=list)
    content: str = ""

about_content = ""
qt_content = ""
download_content = ""

all_pages = [Page("About", "about.html", [Page("Formats", "format.html"),
                                          Page("QtVipster", "gui.html"),
                                          Page("Python API", "python.html")]),
             Page("Download", "download.html"),
            ]

def printPages(pages):
    for page in pages:
        template = env.get_template(page.url)
        with open(str(targetdir.absolute())+"/"+page.url, 'w') as f:
            f.write(template.render(current=page, pages=all_pages, tree=[root]))
        printPages(page.pages)

printPages(all_pages)
