# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

project = 'microkinetics_toolkit'
copyright = '2024, Atsushi Ishikawa'
author = 'Atsushi Ishikawa'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath('..'))

# 必要な拡張機能を読み込む
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',      # 数式(latex)を使う場合
    'sphinx.ext.githubpages',  # githubを使う場合
    'myst_parser',             # markdownを使う場合
]
myst_enable_extensions = [
    "dollarmath",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'  # テーマを変更
html_static_path = ['_static']

html_show_copyright = False
html_show_sphinx = False
