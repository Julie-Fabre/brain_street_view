project = "Brain Street View"
copyright = "2024-2026, Julie M. J. Fabre"
author = "Julie M. J. Fabre"
release = "0.1.1"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = "Brain Street View"
html_static_path = ["_static"]

autodoc_member_order = "bysource"
